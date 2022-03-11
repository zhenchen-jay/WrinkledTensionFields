#pragma once

#include <json/json.h>
#include "WTFSetup.h"
#include "WTFState.h"

namespace WTF
{
	bool loadWTF(const std::string& path, WTFSetup& setup, WTFState& state);
	bool saveWTF(const std::string& path, const WTFSetup& setup, WTFState& state);
}