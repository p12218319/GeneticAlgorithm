/*
	Copyright 2016 Adam Smith - P12218319

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
   email : p12218319@myemail.dmu.ac.uk
*/
#ifndef P12218319_GENETIC_ALGORITHM_BINARY_GENOME_HPP
#define P12218319_GENETIC_ALGORITHM_BINARY_GENOME_HPP

#include <cstdint>

namespace P12218319 { namespace ga{ namespace binary{

    /*!
        \brief
        \detail
        \tparam
        \date 19th January 2016
        \version 0.1
    */
    template<const uint32_t TOTAL_BYTES_>
    struct Genome {
		enum {
			TOTAL_BYTES = TOTAL_BYTES_
		};

		uint8_t genes[TOTAL_BYTES];
    };

}}}

#endif
