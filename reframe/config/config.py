"""ARCHER2 settings"""

from reframe.core.backends import register_launcher
from reframe.core.launchers import JobLauncher

site_configuration = {
    "systems": [
        {
            "name": "archer2",
            "descr": "ARCHER2",
            "hostnames": ["uan", "ln", "dvn"],
            "modules_system": "lmod",
            "partitions": [
                {
                    "name": "login",
                    "descr": "Login nodes",
                    "scheduler": "local",
                    "launcher": "local",
                    "environs": [
                        "Default",
                        "PrgEnv-gnu",
                        "PrgEnv-cray",
                    ],
                },
                {
                    "name": "compute",
                    "descr": "Compute nodes",
                    "scheduler": "slurm",
                    "launcher": "srun",
                    "access": [
                        "--hint=nomultithread",
                        "--distribution=block:block",
                        "--partition=standard",
                        "--qos=standard",
                    ],
                    "environs": [
                        "PrgEnv-gnu",
                        "PrgEnv-cray",
                    ],
                    "max_jobs": 64,
                    "processor": {
                        "num_cpus": 128,
                        "num_cpus_per_socket": 64,
                        "num_sockets": 2,
                    },
                },
                {
                    "name": "compute-gpu",
                    "descr": "Compute nodes with AMD GPUs",
                    "max_jobs": 2,
                    "features": ["gpu"],
                    "scheduler": "slurm",
                    "launcher": "srun",
                    "access": ["--partition=gpu"],
                    "environs": [
                        "rocm-PrgEnv-gnu",
                        "rocm-PrgEnv-cray",
                    ],
                    "resources": [
                        {"name": "qos", "options": ["--qos={qos}"]},
                        {
                            "name": "gpu",
                            "options": ["--gres=gpu:{num_gpus_per_node}"],
                        },
                    ],
                    "processor": {
                        "num_cpus": 32,
                        "num_cpus_per_socket": 32,
                        "num_sockets": 1,
                    },
                    "devices": [
                        {"type": "gpu", "num_devices": 4}
                    ],
                },
            ],
        },
        {
            "name": "cirrus-ex",
            "descr": "Cirrus EX4000",
            "hostnames": ["login"],
            "modules_system": "lmod",
            "partitions": [
                {
                    "name": "login",
                    "descr": "Login nodes",
                    "scheduler": "local",
                    "launcher": "local",
                    "environs": [
                        "Default",
                        "PrgEnv-gnu",
                        "PrgEnv-cray",
                    ],
                },
                {
                    "name": "compute",
                    "descr": "Compute nodes",
                    "scheduler": "slurm",
                    "launcher": "srun",
                    "access": [
                        "--hint=nomultithread",
                        "--distribution=block:block",
                        "--exclusive",
                        "--partition=standard",
                        "--qos=standard",
                    ],
                    "environs": [
                        "PrgEnv-gnu",
                        "PrgEnv-cray",
                    ],
                    "max_jobs": 512,
                    "processor": {
                        "num_cpus": 288,
                        "num_cpus_per_socket": 144,
                        "num_sockets": 2,
                    },
                },
            ],
        }
    ],
    "environments": [
        {
            "name": "PrgEnv-gnu",
            "modules": ["ccs/gnu-2026-06"],
            "cc": "cc",
            "cxx": "CC",
            "ftn": "ftn",
            "target_systems": ["archer2", "cirrus-ex"],
        },
        {
            "name": "PrgEnv-cray",
            "modules": ["ccs/cray-2026-06"],
            "cc": "cc",
            "cxx": "CC",
            "ftn": "ftn",
            "target_systems": ["archer2", "cirrus-ex"],
        },
        {
            "name": "Default",
            "cc": "cc",
            "cxx": "CC",
            "ftn": "ftn",
            "target_systems": ["archer2", "cirrus-ex"],
        },
    ],
    "logging": [
        {
            "level": "debug",
            "perflog_compat": True,
            "handlers": [
                {
                    "type": "stream",
                    "name": "stdout",
                    "level": "info",
                    "format": "%(message)s",
                },
                {
                    "type": "file",
                    "name": "reframe.out",
                    "level": "info",
                    "format": "[%(asctime)s] %(check_info)s: %(message)s",
                    "append": True,
                },
                {
                    "type": "file",
                    "name": "reframe.log",
                    "level": "debug",
                    "format": ("[%(asctime)s] %(levelname)s %(levelno)s: %(check_info)s: %(message)s"),
                    "append": False,
                },
            ],
            "handlers_perflog": [
                {
                    "type": "file",
                    "name": "reframe_perf.out",
                    "level": "info",
                    "perflog_compat": True,
                    "format": (
                        "[%(asctime)s] %(check_info)s: %(check_perf_var)s=%(check_perf_value)s "
                        "(ref=%(check_perf_ref)s;l=%(check_perf_lower_thres)s;u=%(check_perf_upper_thres)s)) "
                        "%(check_perf_unit)s"
                    ),
                    "append": True,
                },
                {
                    "type": "filelog",
                    "prefix": "%(check_system)s/%(check_partition)s",
                    "level": "info",
                    "format": (
                        "%(check_display_name)s|%(check_result)s|%(check_job_completion_time)s|"
                        "%(check_perf_var)s|"
                        "%(check_perf_value)s %(check_perf_unit)s|"
                        "(%(check_perf_ref)s, %(check_perf_lower_thres)s, %(check_perf_upper_thres)s)|"
                    ),
                    "append": True,
                },
            ],
        }
    ],
}
