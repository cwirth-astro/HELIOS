# ==============================================================================
# Module with a database of planetary system parameters used in HELIOS
# Copyright (C) 2020 - 2022 Matej Malik
# ==============================================================================
# This file is part of HELIOS.
#
#     HELIOS is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     HELIOS is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You find a copy of the GNU General Public License in the main
#     HELIOS directory under <license.txt>. If not, see
#     <http://www.gnu.org/licenses/>.
# ==============================================================================

from source import phys_const as pc


class Planet(object):

    def __init__(self, R_p, g_p, a, T_star, R_star, g_star, metal_star, R_p_unit):

        self.R_p = R_p
        if R_p_unit == "R_Earth":
            self.R_p *= pc.R_EARTH / pc.R_JUP

        self.g_p = g_p
        self.a = a
        self.T_star = T_star
        self.R_star = R_star
        self.g_star = g_star
        self.metal_star = metal_star


planet_lib = {}

# Units are [R_p]=R_Earth or R_Jup, [g]=cm s^-2 or [g]=log(cm s^-2), [a]=AU, [R_star]=R_Sun, [T_star]=K

planet_lib["GJ_1214b"] = Planet(R_p=2.85, R_p_unit="R_Earth",
                                g_p=760,
                                a=0.01411,
                                T_star=3026,
                                R_star=0.216,
                                g_star=4.944,
                                metal_star=0.39
                                )  # references: Harpsoe et al. (2013)

planet_lib["HD_209458b"] = Planet(R_p=1.380, R_p_unit="R_Jupiter",
                               g_p=930,
                               a=0.04747,
                               T_star=6117,
                               R_star=1.162,
                               g_star=4.368,
                               metal_star=0.02
                               )  # references: Southworth (2010)

planet_lib["TRAPPIST-1b"] = Planet(R_p=0.0996, R_p_unit="R_Jupiter",
                                g_p=1080,
                                a=0.01154,
                                T_star=2556,
                                R_star=0.121,
                                g_star=5.240,
                                metal_star=0.04
                                )  # references: Agol et al. (2021)

planet_lib["TRAPPIST-1c"] = Planet(R_p=0.0978, R_p_unit="R_Jupiter",
                                g_p=1065,
                                a=0.0158,
                                T_star=2556,
                                R_star=0.121,
                                g_star=5.240,
                                metal_star=0.04
                                )  # references: Agol et al. (2021)

planet_lib["TRAPPIST-1d"] = Planet(R_p=0.0703, R_p_unit="R_Jupiter",
                                g_p=612,
                                a=0.02227,
                                T_star=2556,
                                R_star=0.121,
                                g_star=5.240,
                                metal_star=0.04
                                )  # references: Agol et al. (2021)

planet_lib["TRAPPIST-1e"] = Planet(R_p=0.0821, R_p_unit="R_Jupiter",
                                g_p=801,
                                a=0.02925,
                                T_star=2556,
                                R_star=0.121,
                                g_star=5.240,
                                metal_star=0.04
                                )  # references: Agol et al. (2021)

planet_lib["TRAPPIST-1f"] = Planet(R_p=0.09322, R_p_unit="R_Jupiter",
                                g_p=932,
                                a=0.03849,
                                T_star=2556,
                                R_star=0.121,
                                g_star=5.240,
                                metal_star=0.04
                                )  # references: Agol et al. (2021)

planet_lib["TRAPPIST-1g"] = Planet(R_p=0.1007, R_p_unit="R_Jupiter",
                                g_p=1015,
                                a=0.04683,
                                T_star=2556,
                                R_star=0.121,
                                g_star=5.240,
                                metal_star=0.04
                                )  # references: Agol et al. (2021)

planet_lib["TRAPPIST-1h"] = Planet(R_p=0.06736, R_p_unit="R_Jupiter",
                                g_p=560,
                                a=0.06189,
                                T_star=2556,
                                R_star=0.121,
                                g_star=5.240,
                                metal_star=0.04
                                )  # references: Agol et al. (2021)

planet_lib["WASP-15b"] = Planet(R_p=1.408, R_p_unit="R_Jupiter",
                                g_p=750,
                                a=0.0499,
                                T_star=6372,
                                R_star=1.477,
                                g_star=4.17,
                                metal_star=-0.17
                                )  # references: Bonomo et al. (2017)

planet_lib["KELT-7b"] = Planet(R_p=1.496, R_p_unit="R_Jupiter",
                                g_p=1541,
                                a=0.04415,
                                T_star=6768,
                                R_star=1.81,
                                g_star=4.5,
                                metal_star=0.
                                )  # references: Stassun et al. (2017)

planet_lib["HAT-P-30b"] = Planet(R_p=1.417, R_p_unit="R_Jupiter",
                                g_p=877.4,
                                a=0.0419,
                                T_star=6304,
                                R_star=1.215,
                                g_star=4.36,
                                metal_star=0.13
                                )  # references: Bonomo et al. (2017), Blazek et al. (2022)

planet_lib["NGTS-2b"] = Planet(R_p=1.595, R_p_unit="R_Jupiter",
                                g_p=720,
                                a=0.0630,
                                T_star=6478,
                                R_star=1.702,
                                g_star=4.20,
                                metal_star=-0.06
                                )  # references: Raynard et al. (2018)

planet_lib["TrES-4b"] = Planet(R_p=1.838, R_p_unit="R_Jupiter",
                                g_p=362.5,
                                a=0.05159,
                                T_star=6295,
                                R_star=1.81,
                                g_star=4.09,
                                metal_star=0.28
                                )  # references: Sozzetti et al. (2015)

planet_lib["WASP-94Ab"] = Planet(R_p=1.58, R_p_unit="R_Jupiter",
                                g_p=448.6,
                                a=0.0554,
                                T_star=6170,
                                R_star=1.36,
                                g_star=4.28,
                                metal_star=0.26
                                )  # references: Neveu-VanMalle et al. (2014)

planet_lib["WASP-17b"] = Planet(R_p=1.932, R_p_unit="R_Jupiter",
                                g_p=316.7,
                                a=0.05151,
                                T_star=6650,
                                R_star=1.573,
                                g_star=4.22,
                                metal_star=-0.19
                                )  # references: Anderson et al. (2010)

planet_lib["HD149026b"] = Planet(R_p=0.74, R_p_unit="R_Jupiter",
                                g_p=1720,
                                a=0.0432,
                                T_star=6147,
                                R_star=1.497,
                                g_star=4.20,
                                metal_star=0.36
                                )  # references: Sato et al. (2005)

if __name__ == "__main__":
    print("This module stores information about planetary systems. No guarantee that anything here is remotely correct.")