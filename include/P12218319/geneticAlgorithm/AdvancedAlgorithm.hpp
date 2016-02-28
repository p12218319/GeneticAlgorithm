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
#ifndef P12218319_GENETIC_ALGORITHM_ADVANCED_ALGORITHM_HPP
#define P12218319_GENETIC_ALGORITHM_ADVANCED_ALGORITHM_HPP

#include "Algorithm.hpp"

namespace P12218319 { namespace ga {

    /*!
        \brief
        \detail
        \date 17th January 2016
        \version 0.1
    */
    template<class GENOME, const uint32_t POPULATION_SIZE, const uint32_t CHILDREN_PER_GENERATION, const uint32_t PARENTS_PER_CHILD>
    class P12218319_EXPORT_API AdvancedAlgorithm : public Algorithm<GENOME, POPULATION_SIZE, CHILDREN_PER_GENERATION, PARENTS_PER_CHILD> {
	private:
		typedef Algorithm<GENOME, POPULATION_SIZE, CHILDREN_PER_GENERATION, PARENTS_PER_CHILD> ParentClass;
    private:
        Fitness mMaxFitness;
		Fitness mMinFitness;
		Fitness mAvgFitness;
		uint16_t mMaxAge;
		uint16_t mMinAge;
		uint16_t mAvgAge;
	private:
		void P12218319_CALL PreCalculate() {
			mMinFitness = UINT16_MAX;
            mMaxFitness = 0;
            mAvgFitness = 0;
            mMinAge = UINT16_MAX;
            mMaxAge = 0;
            mAvgAge = 0;
            for(uint32_t i = 0; i < POPULATION_SIZE; ++i) {
                if(mPopulationBuffer[i]->fitness < mMinFitness) mMinFitness = mPopulationBuffer[i]->fitness;
                if(mPopulationBuffer[i]->fitness > mMaxFitness) mMaxFitness = mPopulationBuffer[i]->fitness;
                mAvgFitness += mPopulationBuffer[i]->fitness;

                if(mPopulationBuffer[i]->age < mMinAge) mMinAge = mPopulationBuffer[i]->age;
                if(mPopulationBuffer[i]->age > mMaxAge) mMaxAge = mPopulationBuffer[i]->age;
                mAvgAge += mPopulationBuffer[i]->age;
            }
            mAvgFitness /= POPULATION_SIZE;
            mAvgAge /= POPULATION_SIZE;
		}
    public:
		// Inherited from algorithm
		virtual void P12218319_CALL OnGenerationBegin() override {
			if(ParentClass::GetGeneration() == 0) {
				PreCalculate();
			}
		}

		virtual void P12218319_CALL OnChildGeneration() override {

		}

		virtual void P12218319_CALL OnSurvivorSelection() override {

		}

		virtual void P12218319_CALL OnGenerationEnd() override {
			PreCalculate();
		}

    public:
		P12218319_CALL AdvancedAlgorithm(Randomiser& aRandomiser) :
			Algorithm(aRandomiser)
		{}

        virtual P12218319_CALL ~AdvancedAlgorithm() {

        }

		// population analysis

        inline Fitness P12218319_CALL GetMaxFitness() const {
            return mMaxFitness;
        }

        inline Fitness P12218319_CALL GetMinFitness() const {
            return mMinFitness;
        }

        inline Fitness P12218319_CALL GetAvgFitness() const {
            return mAvgFitness;
        }

        inline uint32_t P12218319_CALL GetMaxAge() const {
            return mMaxAge;
        }

        inline uint32_t P12218319_CALL GetMinAge() const {
            return mMinAge;
        }

        inline uint32_t P12218319_CALL GetAvgAge() const {
            return mAvgAge;
        }

		// Selection helpers

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		void P12218319_CALL SelectRandom(T* aInputs, T* aOutputs) const {
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) {
				aOutputs[i] = aInputs[mRandomiser.Next32u() % INPUT_COUNT];
			}
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		void P12218319_CALL SelectRandomSet(T* aInputs, T* aOutputs) const {
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) {
				bool repeat = true;
				while(repeat) {
					repeat = false;
					aOutputs[i] = aInputs[mRandomiser.Next32u() % INPUT_COUNT];
					for(uint32_t j = 0; j < i; ++j) if (aOutputs[j] == aOutputs[i]) {
						repeat = true;
						break;
					}
				}
			}
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		static void P12218319_CALL SelectFittest(T* aInputs, T* aOutputs) {
			T tmp[INPUT_COUNT];
			for(int32_t i = 0; i < INPUT_COUNT; ++i) tmp[i] = aInputs[i];
			std::sort(tmp, tmp + INPUT_COUNT, [](const ParentClass::Phenotype* const aLeft, const ParentClass::Phenotype* const aRight)->bool {
				return aLeft->fitness > aRight->fitness;
			});
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) aOutputs[i] = tmp[i];
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		static void P12218319_CALL SelectUnfittest(T* aInputs, T* aOutputs) {
			T tmp[INPUT_COUNT];
			for(int32_t i = 0; i < INPUT_COUNT; ++i) tmp[i] = aInputs[i];
			std::sort(tmp, tmp + INPUT_COUNT, [](const ParentClass::Phenotype* const aLeft, const ParentClass::Phenotype* const aRight)->bool {
				return aLeft->fitness < aRight->fitness;
			});
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) aOutputs[i] = tmp[i];
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		static void P12218319_CALL SelectOldest(T* aInputs, T* aOutputs) {
			T tmp[INPUT_COUNT];
			for(int32_t i = 0; i < INPUT_COUNT; ++i) tmp[i] = aInputs[i];
			std::sort(tmp, tmp + INPUT_COUNT, [](const ParentClass::Phenotype* const aLeft, const ParentClass::Phenotype* const aRight)->bool {
				return aLeft->age > aRight->age;
			});
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) aOutputs[i] = tmp[i];
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		static void P12218319_CALL SelectYoungest(T* aInputs, T* aOutputs) {
			T tmp[INPUT_COUNT];
			for(int32_t i = 0; i < INPUT_COUNT; ++i) tmp[i] = aInputs[i];
			std::sort(tmp, tmp + INPUT_COUNT, [](const ParentClass::Phenotype* const aLeft, const ParentClass::Phenotype* const aRight)->bool {
				return aLeft->age < aRight->age;
			});
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) aOutputs[i] = tmp[i];
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		static void P12218319_CALL SelectParents(T* aInputs, T* aOutputs) {
			T tmp[INPUT_COUNT];
			for(int32_t i = 0; i < INPUT_COUNT; ++i) tmp[i] = aInputs[i];
			std::sort(tmp, tmp + INPUT_COUNT, [](const ParentClass::Phenotype* const aLeft, const ParentClass::Phenotype* const aRight)->bool {
				return aLeft->is_parent && ! aRight->is_parent;
			});
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) aOutputs[i] = tmp[i];
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		static void P12218319_CALL SelectNonParents(T* aInputs, T* aOutputs) {
			T tmp[INPUT_COUNT];
			for(int32_t i = 0; i < INPUT_COUNT; ++i) tmp[i] = aInputs[i];
			std::sort(tmp, tmp + INPUT_COUNT, [](const ParentClass::Phenotype* const aLeft, const ParentClass::Phenotype* const aRight)->bool {
				return aRight->is_parent && ! aLeft->is_parent;
			});
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) aOutputs[i] = tmp[i];
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		static void P12218319_CALL SelectChildren(T* aInputs, T* aOutputs) {
			T tmp[INPUT_COUNT];
			for(int32_t i = 0; i < INPUT_COUNT; ++i) tmp[i] = aInputs[i];
			std::sort(tmp, tmp + INPUT_COUNT, [](const ParentClass::Phenotype* const aLeft, const ParentClass::Phenotype* const aRight)->bool {
				return aLeft->is_child && ! aRight->is_child;
			});
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) aOutputs[i] = tmp[i];
		}

		template<const uint32_t INPUT_COUNT, const uint32_t OUTPUT_COUNT, class T = ParentClass::Phenotype*>
		static void P12218319_CALL SelectNonChildren(T* aInputs, T* aOutputs) {
			T tmp[INPUT_COUNT];
			for(int32_t i = 0; i < INPUT_COUNT; ++i) tmp[i] = aInputs[i];
			std::sort(tmp, tmp + INPUT_COUNT, [](const ParentClass::Phenotype* const aLeft, const ParentClass::Phenotype* const aRight)->bool {
				return aRight->is_child && ! aLeft->is_child;
			});
			for(uint32_t i = 0; i < OUTPUT_COUNT; ++i) aOutputs[i] = tmp[i];
		}
    };
}}

#endif
