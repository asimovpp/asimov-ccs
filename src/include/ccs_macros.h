#ifdef VERBOSE
#define dprint(msg) debug_print(msg, __FILE__, __LINE__)
#else
#define dprint(msg) debug_print()
#endif
