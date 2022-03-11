#pragma once
#include "ElasticSetup.h"
#include "ElasticState.h"

bool loadElastic(const std::string& path, ElasticSetup& setup, ElasticState& state);
bool saveElastic(const std::string& path, const ElasticSetup& setup, const ElasticState& state);