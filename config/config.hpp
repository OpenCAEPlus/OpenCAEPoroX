
#ifndef OPENCAEPOROX_CONFIG_HPP
#define OPENCAEPOROX_CONFIG_HPP


#ifdef OCP_CONFIG_FILE
#include OCP_CONFIG_FILE
#else
#ifndef _WIN32
#error "Not found config file."
#endif
#endif



#endif //OPENCAEPOROX_CONFIG_HPP
