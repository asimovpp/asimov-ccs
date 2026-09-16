!v Set of tools to measure flow statistics
module flow_stats

  use mpi

  use core, only: ccs_options
  use kinds, only: ccs_real, ccs_int, CCS_MPI_PRECISION
  use types, only: boundary_patch, ccs_mesh, face_field, field, fluid, cell_locator, &
                   face_locator, neighbour_locator
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi

  use fields, only: get_field
  use meshing, only: count_neighbours, create_cell_locator, create_face_locator, &
                     create_neighbour_locator, get_face_area, get_global_num_cells, &
                     get_boundary_status, get_local_index, get_local_num_cells, get_volume
  use parallel, only: error_handling, is_root
  use timestepping, only: get_current_step, get_current_time, get_timestep, &
                          timestepping_is_active

  implicit none

  private

  !> Boundary-patch metadata augmented with its integrated mass flux.
  type, public, extends(boundary_patch) :: boundary_patch_mass_flux
    real(ccs_real) :: integrated_mass_flux = 0.0_ccs_real
  end type boundary_patch_mass_flux

  !> Integrated mass fluxes and their signed global imbalance.
  type, public :: boundary_mass_flux_balance
    type(boundary_patch_mass_flux), dimension(:), allocatable :: patches
    real(ccs_real) :: global_imbalance = 0.0_ccs_real
  end type boundary_mass_flux_balance

  !> Boundary metadata and owned face geometry reused across mass-flux integrations.
  !> Rebuild after changes to boundary metadata, mesh topology, geometry, or face-field indexing.
  type, public :: boundary_flux_context
    private
    type(boundary_patch), dimension(:), allocatable :: patch_metadata
    integer(ccs_int), dimension(:), allocatable :: patch_ids
    integer(ccs_int), dimension(:), allocatable :: indices_f
    real(ccs_real), dimension(:), allocatable :: areas_f
  end type boundary_flux_context

  public :: create_boundary_flux_context
  public :: integrate_boundary_mass_fluxes
  public :: report_boundary_mass_flux_balance
  public :: report_cfl

contains

  !> Collect the metadata and geometry needed to integrate mass flux through each boundary patch.
  !> The supplied mesh must be the active mesh used by the meshing accessors. Only faces of
  !> locally owned cells are stored, so each physical boundary face contributes on one rank.
  !> The returned context stores patch IDs, face-field indices, and areas, but no flux values
  !> or MPI state. Rebuild it if boundary metadata, topology, geometry, or face indexing changes;
  !> changes to the mass-flux values alone do not require rebuilding it.
  pure function create_boundary_flux_context(mesh) result(context)
    type(ccs_mesh), intent(in) :: mesh
    type(boundary_flux_context) :: context

    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    integer(ccs_int) :: index_p
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: index_f
    integer(ccs_int) :: nnb
    integer(ccs_int) :: i
    integer(ccs_int) :: k
    integer(ccs_int) :: patch_id
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: boundary_face_count
    logical :: is_boundary

    if (.not. allocated(mesh%boundary_patches)) then
      error stop "Boundary patch metadata is not allocated"
    end if
    call get_local_num_cells(local_num_cells)
    if (local_num_cells < 0_ccs_int) then
      error stop "Local cell count must be nonnegative"
    end if

    ! Preserve boundary-ID order so all ranks accumulate into corresponding patch entries.
    context%patch_metadata = mesh%boundary_patches
    ! Count owned physical boundary faces before allocating the cache to its exact size.
    boundary_face_count = 0_ccs_int
    do index_p = 1, local_num_cells
      call create_cell_locator(index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      do i = 1, nnb
        call create_neighbour_locator(loc_p, i, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        if (is_boundary) boundary_face_count = boundary_face_count + 1_ccs_int
      end do
    end do

    allocate (context%patch_ids(boundary_face_count))
    allocate (context%indices_f(boundary_face_count))
    allocate (context%areas_f(boundary_face_count))

    ! Visit the same faces again to store the geometry needed by subsequent integrations.
    k = 0_ccs_int
    do index_p = 1, local_num_cells
      call create_cell_locator(index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      do i = 1, nnb
        call create_neighbour_locator(loc_p, i, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        if (is_boundary) then
          call get_local_index(loc_nb, index_nb)
          ! Physical boundary neighbours encode the boundary ID as a negative index.
          patch_id = -index_nb
          if (patch_id < 1_ccs_int .or. &
              patch_id > size(context%patch_metadata, kind=ccs_int)) then
            error stop "Boundary face has an invalid patch ID"
          end if
          k = k + 1_ccs_int
          call create_face_locator(index_p, i, loc_f)
          call get_local_index(loc_f, index_f)
          ! Keep the field index and area together; the current flux is read at integration time.
          context%patch_ids(k) = patch_id
          context%indices_f(k) = index_f
          call get_face_area(loc_f, context%areas_f(k))
        end if
      end do
    end do
  end function create_boundary_flux_context

  !> Integrate the current mass flux over each boundary patch and calculate the global imbalance.
  !> Each owned boundary face contributes its mass-flux value multiplied by its area. Contributions
  !> are summed across the communicator, giving every rank identical totals in boundary-ID order.
  !> Positive flux is outward and negative flux is inward; the imbalance is the signed sum of patches.
  !> All ranks must call this function with matching patch metadata. The context must still describe
  !> the mesh geometry and face-field indexing for mass_flux_f (see create_boundary_flux_context).
  function integrate_boundary_mass_fluxes(par_env, mass_flux_f, context) result(balance)
    class(parallel_environment), intent(in) :: par_env
    type(face_field), intent(in) :: mass_flux_f
    type(boundary_flux_context), intent(in) :: context
    type(boundary_mass_flux_balance) :: balance

    real(ccs_real), dimension(:), allocatable :: global_patch_fluxes
    real(ccs_real), dimension(:), allocatable :: local_patch_fluxes
    integer :: ierr
    integer :: patch_count

    call validate_boundary_flux_inputs(mass_flux_f, context)

    ! Accumulate local contributions, including zero entries for patches with no owned faces.
    local_patch_fluxes = integrate_local_boundary_fluxes( &
                         size(context%patch_metadata, kind=ccs_int), &
                         context%patch_ids, context%indices_f, context%areas_f, &
                         mass_flux_f%values_ro)

    patch_count = size(context%patch_metadata)
    allocate (global_patch_fluxes(patch_count), source=0.0_ccs_real)
    ! Sum corresponding patch entries across all ranks and return every total to every rank.
    ! Ranks without boundary faces contribute zeros; an empty patch list needs no collective.
    if (patch_count > 0) then
      select type (par_env)
      type is (parallel_environment_mpi)
        call MPI_Allreduce(local_patch_fluxes, global_patch_fluxes, patch_count, &
                           CCS_MPI_PRECISION, MPI_SUM, par_env%comm, ierr)
        call error_handling(ierr, "mpi", par_env)
      class default
        error stop "Unsupported parallel environment"
      end select
    end if

    balance = assemble_boundary_flux_balance(context%patch_metadata, global_patch_fluxes)
  end function integrate_boundary_mass_fluxes

  !> Check that the context and current face-field storage can be indexed safely during integration.
  !> Require allocated caches with matching sizes, valid patch IDs, associated read-only flux values,
  !> and face indices within their extent. These checks do not detect changes to cached mesh geometry.
  pure subroutine validate_boundary_flux_inputs(mass_flux_f, context)
    type(face_field), intent(in) :: mass_flux_f
    type(boundary_flux_context), intent(in) :: context

    if (.not. allocated(context%patch_metadata)) then
      error stop "Boundary flux context is not initialized"
    end if
    if (.not. allocated(context%patch_ids) .or. &
        .not. allocated(context%indices_f) .or. &
        .not. allocated(context%areas_f)) then
      error stop "Boundary flux context cache is not allocated"
    end if
    if (size(context%patch_ids) /= size(context%indices_f) .or. &
        size(context%patch_ids) /= size(context%areas_f)) then
      error stop "Boundary flux context cache extents do not match"
    end if
    if (size(context%patch_ids) > 0) then
      if (minval(context%patch_ids) < 1_ccs_int .or. &
          maxval(context%patch_ids) > size(context%patch_metadata, kind=ccs_int)) then
        error stop "Cached boundary patch ID is invalid"
      end if
    end if

    if (.not. associated(mass_flux_f%values_ro)) then
      error stop "Boundary mass-flux read-only values are not associated"
    end if
    if (size(context%indices_f) > 0) then
      if (minval(context%indices_f) < 1_ccs_int .or. &
          maxval(context%indices_f) > size(mass_flux_f%values_ro, kind=ccs_int)) then
        error stop "Cached boundary face index exceeds the face-field extent"
      end if
    end if
  end subroutine validate_boundary_flux_inputs

  !> Sum mass-flux value times face area for each patch using this rank's owned boundary faces.
  !> Patches without local faces retain zero; no communication or sign adjustment is performed here.
  pure function integrate_local_boundary_fluxes(patch_count, patch_ids, indices_f, &
                                                areas_f, mass_flux_values) result(patch_fluxes)
    integer(ccs_int), intent(in) :: patch_count
    integer(ccs_int), dimension(:), intent(in) :: patch_ids
    integer(ccs_int), dimension(:), intent(in) :: indices_f
    real(ccs_real), dimension(:), intent(in) :: areas_f
    real(ccs_real), dimension(:), intent(in) :: mass_flux_values
    real(ccs_real), dimension(patch_count) :: patch_fluxes

    integer(ccs_int) :: i

    patch_fluxes = 0.0_ccs_real
    do i = 1, size(patch_ids, kind=ccs_int)
      patch_fluxes(patch_ids(i)) = patch_fluxes(patch_ids(i)) + &
                                   mass_flux_values(indices_f(i)) * areas_f(i)
    end do
  end function integrate_local_boundary_fluxes

  !> Combine boundary metadata and globally reduced fluxes in boundary-ID order.
  !> The global imbalance is their signed sum; an empty patch list gives zero imbalance.
  pure function assemble_boundary_flux_balance(patch_metadata, global_patch_fluxes) result(balance)
    type(boundary_patch), dimension(:), intent(in) :: patch_metadata
    real(ccs_real), dimension(:), intent(in) :: global_patch_fluxes
    type(boundary_mass_flux_balance) :: balance

    integer(ccs_int) :: i

    allocate (balance%patches(size(patch_metadata)))
    do i = 1, size(patch_metadata, kind=ccs_int)
      balance%patches(i)%name = patch_metadata(i)%name
      balance%patches(i)%patch_type = patch_metadata(i)%patch_type
      balance%patches(i)%integrated_mass_flux = global_patch_fluxes(i)
    end do
    balance%global_imbalance = sum(global_patch_fluxes)
  end function assemble_boundary_flux_balance

  !> Print an already integrated balance on communicator rank zero when diagnostics output is enabled.
  !> Use a step/time heading for an active transient solve and a steady-solve heading otherwise.
  !> The diagnostics flag controls printing only: the solver integrates the balance every step.
  subroutine report_boundary_mass_flux_balance(par_env, run_options, balance)
    use logging, only: log_unit_out

    class(parallel_environment), intent(in) :: par_env
    type(ccs_options), intent(in) :: run_options
    type(boundary_mass_flux_balance), intent(in) :: balance

    character(len=14) :: time_text
    integer(ccs_int) :: step
    real(ccs_real) :: time

    if (.not. run_options%diagnostics%boundary_mass_fluxes) return
    if (par_env%proc_id /= 0_ccs_int) return

    if (timestepping_is_active()) then
      call get_current_step(step)
      call get_current_time(time)
      write (time_text, '(es14.6)') time
      write (log_unit_out, '("Boundary mass-flux balance: step ",i0,", time ",a)') &
        step, trim(adjustl(time_text))
    else
      write (log_unit_out, '(a)') "Boundary mass-flux balance: steady solve"
    end if
    write (log_unit_out, '(a)') "Sign convention: positive outward, negative inward"
    call write_boundary_mass_flux_table(balance)
  end subroutine report_boundary_mass_flux_balance

  !> Write an aligned table of signed patch fluxes and their global imbalance to the log output.
  !> Size columns from the headers and patch metadata, and identify an empty patch list explicitly.
  !> The caller is responsible for restricting output to the reporting rank.
  subroutine write_boundary_mass_flux_table(balance)
    use logging, only: log_unit_out

    type(boundary_mass_flux_balance), intent(in) :: balance

    character(len=*), parameter :: id_header = "Boundary ID"
    character(len=*), parameter :: name_header = "Boundary name"
    character(len=*), parameter :: type_header = "Patch type"
    character(len=*), parameter :: flux_header = "Mass flux [solver units]"
    character(len=*), parameter :: global_label = "GLOBAL IMBALANCE"
    integer :: flux_width
    integer :: i
    integer :: id_width
    integer :: name_width
    integer :: name_column
    integer :: numeric_flux_column
    integer :: separator_one_column
    integer :: separator_two_column
    integer :: separator_three_column
    integer :: separator_width
    integer :: type_column
    integer :: type_width
    integer :: flux_column
    character(len=32) :: id_text
    character(len=256) :: global_format
    character(len=256) :: header_format
    character(len=256) :: patch_format
    character(len=:), allocatable :: separator

    ! Allow both the headers and all patch names/types to fit without truncation.
    id_width = len(id_header)
    name_width = max(len(name_header), len(global_label))
    type_width = len(type_header)
    flux_width = max(len(flux_header), 14)
    do i = 1, size(balance%patches)
      write (id_text, '(i0)') i
      id_width = max(id_width, len_trim(id_text))
      name_width = max(name_width, len_trim(balance%patches(i)%name))
      type_width = max(type_width, len_trim(balance%patches(i)%patch_type))
    end do
    separator_width = id_width + name_width + type_width + flux_width + 9
    separator = repeat("-", separator_width)

    ! Convert the column widths to positions shared by the header, patch, and imbalance rows.
    separator_one_column = id_width + 1
    name_column = separator_one_column + 3
    separator_two_column = name_column + name_width
    type_column = separator_two_column + 3
    separator_three_column = type_column + type_width
    flux_column = separator_three_column + 3
    numeric_flux_column = flux_column + flux_width - 14

    ! Build formats for the calculated column positions, retaining explicit signs on fluxes.
    write (header_format, &
           '("(A,T",I0,",A,T",I0,",A,T",I0,",A,T",I0,",A,T",I0,",A,T",I0,",A)")') &
      separator_one_column, name_column, separator_two_column, type_column, &
      separator_three_column, flux_column
    write (patch_format, &
           '("(I0,T",I0,",A,T",I0,",A,T",I0,",A,T",I0,",A,T",I0,",A,T",I0,",SP,ES14.6)")') &
      separator_one_column, name_column, separator_two_column, type_column, &
      separator_three_column, numeric_flux_column
    write (global_format, &
           '("(T",I0,",A,T",I0,",A,T",I0,",A,T",I0,",A,T",I0,",SP,ES14.6)")') &
      separator_one_column, name_column, separator_two_column, separator_three_column, &
      numeric_flux_column

    write (log_unit_out, '(a)') separator
    write (log_unit_out, trim(header_format)) id_header, " | ", name_header, " | ", &
      type_header, " | ", flux_header
    write (log_unit_out, '(a)') separator
    if (size(balance%patches) == 0) then
      write (log_unit_out, '(a)') "(no boundary patches)"
    else
      do i = 1, size(balance%patches)
        write (log_unit_out, trim(patch_format)) i, " | ", &
          trim(balance%patches(i)%name), " | ", &
          trim(balance%patches(i)%patch_type), " | ", &
          balance%patches(i)%integrated_mass_flux
      end do
    end if
    write (log_unit_out, '(a)') separator
    write (log_unit_out, trim(global_format)) " | ", global_label, " | ", " | ", &
      balance%global_imbalance
    write (log_unit_out, '(a)') separator
  end subroutine write_boundary_mass_flux_table

  subroutine report_cfl(par_env, flow)
    use logging, only: log_unit_out
    class(parallel_environment), intent(in) :: par_env
    type(fluid), intent(in) :: flow

    class(field), pointer :: u, v, w

    real(ccs_real) :: cfl_max, cfl_avg, cfl_i
    real(ccs_real) :: dt
    real(ccs_real) :: dx, V_p
    type(cell_locator) :: loc_p

    real(ccs_real) :: vel

    integer(ccs_int) :: i, nlocal, nglobal

    integer :: ierr

    if (.not. timestepping_is_active()) then
      return
    end if

    dt = get_timestep()

    cfl_max = 0.0_ccs_real
    cfl_avg = 0.0_ccs_real

    call get_field(flow, "u", u)
    call get_field(flow, "v", v)
    call get_field(flow, "w", w)

    call get_local_num_cells(nlocal)
    call get_global_num_cells(nglobal)
    do i = 1, nlocal
      vel = norm2([u%values_ro(i), v%values_ro(i), w%values_ro(i)])

      call create_cell_locator(i, loc_p)
      call get_volume(loc_p, V_p)
      dx = get_lscale(V_p)

      cfl_i = cfl(vel, dt, dx)
      cfl_max = max(cfl_i, cfl_max)
      cfl_avg = cfl_avg + cfl_i / nglobal
    end do

    nullify (u, v, w)

    ! Reduce the MAX/AVG CFL numbers
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, cfl_max, 1, MPI_DOUBLE, MPI_MAX, par_env%comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, cfl_avg, 1, MPI_DOUBLE, MPI_SUM, par_env%comm, ierr)
    class default
      error stop "Unsupported parallel environment"
    end select

    if (is_root(par_env)) then
      write (log_unit_out, *) "CFL Max: ", cfl_max
      write (log_unit_out, *) "CFL Avg: ", cfl_avg
    end if

  end subroutine report_cfl

  pure real(ccs_real) function get_lscale(vol) result(l)
    real(ccs_real), intent(in) :: vol

    l = vol**(1.0 / 3.0)
  end function get_lscale

  pure real(ccs_real) function cfl(u, dt, dx)
    real(ccs_real), intent(in) :: u
    real(ccs_real), intent(in) :: dt
    real(ccs_real), intent(in) :: dx

    cfl = u * dt / dx

  end function cfl

end module flow_stats
