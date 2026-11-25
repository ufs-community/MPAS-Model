#include "ezxml_shim.h"

ezxml_t shim_parse_file(const char *file) {
    return ezxml_parse_file(file);
}

ezxml_t shim_child(ezxml_t xml, const char *name) {
    return ezxml_child(xml, name);
}

ezxml_t shim_next(ezxml_t xml) {
    return xml ? xml->next : (ezxml_t)0;
}

const char *shim_name(ezxml_t xml) {
    return xml ? xml->name : (const char *)0;
}

const char *shim_txt(ezxml_t xml) {
    return xml ? xml->txt : "";
}

const char *shim_attr(ezxml_t xml, const char *attrname) {
    return ezxml_attr(xml, attrname);
}

void shim_free(ezxml_t xml) {
    ezxml_free(xml);
}
