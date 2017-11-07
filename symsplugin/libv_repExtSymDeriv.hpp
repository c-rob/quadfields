
#pragma once

#ifdef _WIN32
    #define VREP_DLLEXPORT extern "C" __declspec(dllexport)
#else
    #define VREP_DLLEXPORT extern "C"
#endif


#include "v_repLib.h"
#include <iostream>

#ifdef _WIN32
    #ifdef QT_COMPIL
        #include <direct.h>
    #else
        #include <shlwapi.h>
        #pragma comment(lib, "Shlwapi.lib")
    #endif
#else
    #include <unistd.h>
#endif

#ifdef __APPLE__
#define _stricmp strcmp
#endif

#define PLUGIN_VERSION 4 // 2 since version 3.2.1, 3 since V3.3.1, 4 since V3.4.0


// The 3 required entry points of the V-REP plugin:
VREP_DLLEXPORT unsigned char v_repStart(void* reservedPointer,int reservedInt);
VREP_DLLEXPORT void v_repEnd();
VREP_DLLEXPORT void* v_repMessage(int message,int* auxiliaryData,void* customData,int* replyData);
