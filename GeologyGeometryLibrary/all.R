


# Copyright 2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# This file is nothing more than a convenient way to load our entire R library and its dependencies (excluding the C part).

# Load external libraries. For example, expm must be loaded before def.R, because the latter defines defExp to be expm.
library("rgl")
library("fields")
library("MASS")
library("ICSNP")
library("expm")
library("Directional")
library("pracma")
library("ggplot2")

# Load our geologyGeometry library, excluding the C part.
source("GeologyGeometryLibrary/miscellany.R")
source("GeologyGeometryLibrary/ray.R")
source("GeologyGeometryLibrary/rayFisher.R")
source("GeologyGeometryLibrary/rayFisherVw.R")
source("GeologyGeometryLibrary/rayWellner.R")
source("GeologyGeometryLibrary/rayRegression.R")
source("GeologyGeometryLibrary/rayPlot.R")
source("GeologyGeometryLibrary/ray2D.R")
source("GeologyGeometryLibrary/line.R")
source("GeologyGeometryLibrary/lineUniform.R")
source("GeologyGeometryLibrary/lineWatson.R")
source("GeologyGeometryLibrary/lineBingham.R")
source("GeologyGeometryLibrary/lineRegression.R")
source("GeologyGeometryLibrary/lineWellner.R")
source("GeologyGeometryLibrary/linePlot.R")
source("GeologyGeometryLibrary/line2D.R")
source("GeologyGeometryLibrary/rot.R")
source("GeologyGeometryLibrary/rotRancourt.R")
source("GeologyGeometryLibrary/rotUniform.R")
source("GeologyGeometryLibrary/rotFisher.R")
source("GeologyGeometryLibrary/rotRegression.R")
source("GeologyGeometryLibrary/rotPlot.R")
source("GeologyGeometryLibrary/ori.R")
source("GeologyGeometryLibrary/oriRegression.R")
source("GeologyGeometryLibrary/oriPlot.R")
source("GeologyGeometryLibrary/ell.R")
source("GeologyGeometryLibrary/ellSPO.R")
source("GeologyGeometryLibrary/ellPlot.R")
source("GeologyGeometryLibrary/def.R")
source("GeologyGeometryLibrary/defHomogeneous.R")
source("GeologyGeometryLibrary/defJeffery.R")
source("GeologyGeometryLibrary/defEshelby.R")
source("GeologyGeometryLibrary/geology.R")

