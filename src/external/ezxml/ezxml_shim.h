#ifndef EZXML_SHIM_H
#define EZXML_SHIM_H

#include "ezxml.h"

#ifdef __cplusplus
extern "C" {
#endif

ezxml_t     shim_parse_file(const char *file);
ezxml_t     shim_child(ezxml_t xml, const char *name);
ezxml_t     shim_next(ezxml_t xml);       /* next element with the SAME tag name */
const char *shim_name(ezxml_t xml);
const char *shim_txt(ezxml_t xml);
const char *shim_attr(ezxml_t xml, const char *attrname);
void        shim_free(ezxml_t xml);

#ifdef __cplusplus
}
#endif

#endif /* EZXML_SHIM_H */
