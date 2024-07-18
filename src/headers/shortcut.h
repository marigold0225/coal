//
// Created by mafu on 4/2/2024.
//

#pragma once
#include "particle.h"
#include <map>
#include <vector>

namespace Coal {

    using Pair               = std::pair<int, int>;
    using RapidityRange      = std::pair<double, double>;
    using RapidityArray      = std::vector<RapidityRange>;
    using ParticleArray      = std::vector<Particle>;
    using MultiParticleArray = std::vector<ParticleArray>;
    //[pdg,particles]
    using ParticleTypeMap = std::map<int, ParticleArray>;
    //[eventID,particles]
    using ParticleEventMap = std::map<int, ParticleArray>;
    using EventsMap        = std::map<int, ParticleTypeMap>;
    //[centrility, EventMap]
    using CentralityMap = std::map<Pair, EventsMap>;
    //[eventID,EventPlane angle]
    using EventPlaneMap      = std::map<int, double>;
    using FlowType = std::map<std::string, double>;
    using EventPlaneMapArray = std::map<int, std::map<std::string,double>>;

}// namespace Coal
