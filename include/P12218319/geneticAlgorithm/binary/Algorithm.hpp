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
#ifndef P12218319_GENETIC_ALGORITM_BINARY_ALGORITHM_HPP
#define P12218319_GENETIC_ALGORITM_BINARY_ALGORITHM_HPP

#include "..\AdvancedAlgorithm.hpp"

namespace P12218319 { namespace ga { namespace binary{

    /*!
        \brief
        \detail
        \date 19th January 2016
        \version 0.1
    */
    template<class GENOME, const uint32_t POPULATION_SIZE, const uint32_t CHILDREN_PER_GENERATION, const uint32_t PARENTS_PER_CHILD>
    class P12218319_EXPORT_API Algorithm : public AdvancedAlgorithm<GENOME, POPULATION_SIZE, CHILDREN_PER_GENERATION, PARENTS_PER_CHILD> {
	public:
		enum MutationMode {
			MODE_SET,
			MODE_CLEAR,
			MODE_FLIP,
			MODE_RANDOM
		};
    private:
        typedef AdvancedAlgorithm<GENOME, POPULATION_SIZE, CHILDREN_PER_GENERATION, PARENTS_PER_CHILD> ParentClass;
    public:
		P12218319_CALL Algorithm(Randomiser& aRandomiser) :
			AdvancedAlgorithm(aRandomiser)
		{}

        virtual P12218319_CALL ~Algorithm() {

        }

		// mutation

		template<const uint8_t MUTATION_RATE, const MutationMode MODE = MODE_FLIP, const uint32_t MAX_MUTATIONS = ParentClass::Genome::TOTAL_BYTES * 8>
		void P12218319_CALL RandomMutationBit(typename ParentClass::Genome& aGenome) const {
			// Perform mutation
			uint32_t mutations = 0;
			for (uint32_t i = 0; i < ParentClass::Genome::TOTAL_BYTES; ++i) {
				uint8_t mask = 1;
				for(uint32_t j = 0; j < 8; ++j) {
					if((mRandomiser.Next32u() % 100) < MUTATION_RATE) {
						switch(MODE) {
						case MODE_SET:
							aGenome.genes[i] |= mask;
							break;
						case MODE_CLEAR:
							aGenome.genes[i] &= ~mask;
							break;
						case MODE_FLIP:
							if(aGenome.genes[i] & mask) {
								aGenome.genes[i] &= ~mask;
							}else {
								aGenome.genes[i] |= ~mask;
							}
							break;
						case MODE_RANDOM:
							aGenome.genes[i] &= ~mask;
							aGenome.genes[i] |= (mRandomiser.Next32u() % 100) > 50 ? 1 : 0;
							break;
						}
						++mutations;
						if(mutations == MAX_MUTATIONS) return;
					}
					mask <<= 1;
				}
			}
		}

		template<const uint8_t MUTATION_RATE, const MutationMode MODE = MODE_RANDOM, const uint32_t MAX_MUTATIONS = ParentClass::Genome::TOTAL_BYTES>
		void P12218319_CALL RandomMutationByte(typename ParentClass::Genome& aGenome) const {
			uint32_t mutations = 0;
			for (uint32_t i = 0; i < ParentClass::Genome::TOTAL_BYTES; ++i) {
				if ((mRandomiser.Next32u() % 100) < MUTATION_RATE) {
					switch(MODE) {
					case MODE_SET:
						aGenome.genes[i] = UINT8_MAX;
						break;
					case MODE_CLEAR:
						aGenome.genes[i] = 0;
						break;
					case MODE_FLIP:
						aGenome.genes[i] = ~ aGenome.genes[i];
						break;
					case MODE_RANDOM:
						aGenome.genes[i] = mRandomiser.Next32u() % UINT8_MAX;
						break;
					}
				}
				++mutations;
				if (mutations == MAX_MUTATIONS) return;
			}
		}

		template<const uint8_t MUTATION_RATE, const uint32_t SWAPS = 1>
		void P12218319_CALL SwapBits(typename ParentClass::Genome& aGenome) const {
			for(uint32_t i = 0; i < SWAPS; ++i) if((mRandomiser.Next32u() % 100) < MUTATION_RATE) {
				const uint32_t bit0 = mRandomiser.Next32u() % (ParentClass::Genome::TOTAL_BYTES * 8);
				const uint32_t bit1 = mRandomiser.Next32u() % (ParentClass::Genome::TOTAL_BYTES * 8);
				const uint32_t byte0 = bit0 / 8;
				const uint32_t byte1 = bit1 / 8;
				const uint32_t offset0 = 1 << (byte0 % 8);
				const uint32_t offset1 = 1 << (byte1 % 8);

				const uint8_t value0 = aGenome.genes[byte0] & offset0;
				const uint8_t value1 = aGenome.genes[byte1] & offset1;
				aGenome.genes[byte0] &= offset0;
				aGenome.genes[byte1] &= offset1;
				aGenome.genes[byte0] |= value1 ? offset0 : 0;
				aGenome.genes[byte1] |= value0 ? offset1 : 0;
			}
		}

		// crossover

		void P12218319_CALL UniformCrossoverBit(typename ParentClass::Genome& aGenome) const {
			for(uint32_t i = 0; i < ParentClass::Genome::TOTAL_BYTES; ++i) {
				uint8_t mask = 1;
				aGenome.genes[i] = 0;
				for(uint32_t j = 0; j < 8; ++j) {
					const uint8_t* const parentGenome = mParentBuffer[mRandomiser.Next32u() % ParentClass::PARENTS_PER_CHILD]->genome.genes;
					aGenome.genes[i] |= parentGenome[i] & mask;
					mask <<= 1;
				}
			}
		}

		void P12218319_CALL UniformCrossoverByte(typename ParentClass::Genome& aGenome) const {
			for(uint32_t i = 0; i < ParentClass::Genome::TOTAL_BYTES; ++i) {
				aGenome.genes[i] = mParentBuffer[mRandomiser.Next32u() % ParentClass::PARENTS_PER_CHILD]->genome.genes[i];
			}
		}

		void P12218319_CALL HighByteCrossover(typename ParentClass::Genome& aGenome) const {
			for(uint32_t i = 0; i < ParentClass::Genome::TOTAL_BYTES; ++i) {
				uint8_t high = mParentBuffer[0]->genes[i];
				for(uint32_t j = 1; j < ParentClass::PARENTS_PER_CHILD; ++j) if(high < mParentBuffer[j]->genes[i]) high = mParentBuffer[j]->genome.genes[i];
				aGenome.genes[i] = high;
			}
		}

		void P12218319_CALL LowByteCrossover(typename ParentClass::Genome& aGenome) const {
			for (uint32_t i = 0; i < ParentClass::Genome::TOTAL_BYTES; ++i) {
				uint8_t low = mParentBuffer[0]->genes[i];
				for(uint32_t j = 1; j < ParentClass::PARENTS_PER_CHILD; ++j) if (low > mParentBuffer[j]->genes[i]) low = mParentBuffer[j]->genome.genes[i];
				aGenome.genes[i] = low;
			}
		}

		void P12218319_CALL OnePointCrossoverByte(typename ParentClass::Genome& aGenome) const {
			enum {
				avg_length = ParentClass::Genome::TOTAL_BYTES / ParentClass::PARENTS_PER_CHILD
			};

			// Calculate crossover points
			uint16_t crossoverLengths[ParentClass::PARENTS_PER_CHILD];
			{
				uint32_t total = 0;
				for(uint32_t i = 0; i < ParentClass::PARENTS_PER_CHILD; ++i) {
					crossoverLengths[i] = mRandomiser.Next32u() % avg_length;
					total += crossoverLengths[i];
				}
				const uint32_t leftover = ParentClass::Genome::TOTAL_BYTES - total;
				const uint32_t add = leftover / ParentClass::PARENTS_PER_CHILD;
				for(uint32_t i = 0; i < ParentClass::PARENTS_PER_CHILD; ++i) crossoverLengths[i] += add;
			}

			// Perform crossover
			uint32_t offset = 0;
			for(uint32_t i = 0; i < ParentClass::PARENTS_PER_CHILD; ++i) {
				for(uint32_t j = 0; j < crossoverLengths[i]; ++j) aGenome.genes[offset + j] = mParentBuffer[i]->genome.genes[offset + j];
					offset += crossoverLengths[i];
			}
		}
    };
}}}

#endif
