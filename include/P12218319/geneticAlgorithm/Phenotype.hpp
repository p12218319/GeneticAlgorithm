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
#ifndef P12218319_GENETIC_ALGORITHM_PHENOTYPE_HPP
#define P12218319_GENETIC_ALGORITHM_PHENOTYPE_HPP

#include <cstdint>

namespace P12218319 { namespace ga{

	typedef uint16_t Fitness;

    /*!
        \brief
        \detail
        \tparam GENOME The genome type for this individual
        \date 17th January 2016
        \version 0.1
    */
	template<class GENOME>
	struct Phenotype {
		typedef GENOME Genome;     //!< The genome type for this phenotype
		
		Genome genome;             //!< The genome that contains the genetic information
		struct {
			uint8_t isParent   : 1;    //!< Set to 1 when this phenotype has 1 or more children in the current generation
			uint8_t isChild    : 1;    //!< Set to 1 if this phenotype was created as a child during the current generation
		};
		uint16_t age;                   //!< The number of generations the phenotype has survived
		Fitness fitness;				//!< Stores the fitness of the phenotype when determined
    };

}}

#endif
