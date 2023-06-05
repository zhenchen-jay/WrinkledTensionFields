#pragma once

#include "TFWSetup.h"
#include "TFWState.h"

namespace TFW
{
	bool loadTFW(const std::string& path, TFWSetup& setup, TFWState& state);
	bool saveTFW(const std::string& path, const TFWSetup& setup, TFWState& state);
}