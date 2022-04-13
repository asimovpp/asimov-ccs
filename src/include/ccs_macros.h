#ifdef VERBOSE
#define diag(msg) diag_print(msg, __FILE__, __LINE__)
#else
#define diag(msg) diag_null(msg)
#endif
