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
#ifndef P12218319_GENETIC_ALGORITHM_ALGORITHM_HPP
#define P12218319_GENETIC_ALGORITHM_ALGORITHM_HPP

#include "Phenotype.hpp"
#include "P12218319\core\Randomiser.hpp"

namespace P12218319 { namespace ga {

    /*!
        \brief
        \detail
        \date 17th January 2016
        \version 0.1
    */
    template<class GENOME, const uint32_t POPULATION_SIZE_, const uint32_t CHILDREN_PER_GENERATION_, const uint32_t PARENTS_PER_CHILD_>
    class P12218319_EXPORT_API Algorithm {
    public:
        enum {
            POPULATION_SIZE = POPULATION_SIZE_,
            CHILDREN_PER_GENERATION = CHILDREN_PER_GENERATION_,
            PARENTS_PER_CHILD = PARENTS_PER_CHILD_
        };
        typedef GENOME Genome;
        typedef Phenotype<Genome> Phenotype;
    private:
		Phenotype mPopulationBufferA[POPULATION_SIZE];
		Phenotype mPopulationBufferB[POPULATION_SIZE];
        uint32_t mGenerationCount;
	protected:
		Randomiser& mRandomiser;
		Phenotype** mPopulationBuffer;
		Phenotype** mChildBuffer;
		Phenotype** mParentBuffer;
		Phenotype** mSurvivorBuffer;
    protected:
        virtual void P12218319_CALL Initialise(Genome& aGenome) const = 0;
        virtual Fitness P12218319_CALL CalculateFitness(const Genome& aGenome) const = 0;
        virtual void P12218319_CALL SelectParents() = 0;
        virtual void P12218319_CALL Crossover(Genome& aGenome) const = 0;
        virtual void P12218319_CALL Mutate(Genome& aGenome) const = 0;
        virtual void P12218319_CALL SelectSurvivors() = 0;
        virtual bool P12218319_CALL ShouldTerminate() const = 0;

		virtual void P12218319_CALL OnGenerationBegin() = 0;
		virtual void P12218319_CALL OnChildGeneration() = 0;
		virtual void P12218319_CALL OnSurvivorSelection() = 0;
		virtual void P12218319_CALL OnGenerationEnd() = 0;
    public:
		P12218319_CALL Algorithm(Randomiser& aRandomiser) :
            mGenerationCount(0),
			mRandomiser(aRandomiser)
        {}

        virtual P12218319_CALL ~Algorithm() {
        }

		inline uint32_t P12218319_CALL GetGeneration() const {
			return mGenerationCount;
		}

        const Phenotype* P12218319_CALL operator()() {
			Phenotype* population = mPopulationBufferA;
			Phenotype children[CHILDREN_PER_GENERATION];
			Phenotype* populationBuffer[POPULATION_SIZE + CHILDREN_PER_GENERATION];
			Phenotype* parentBuffer[PARENTS_PER_CHILD];
			Phenotype* survivorBuffer[POPULATION_SIZE];

			mPopulationBuffer = populationBuffer;
			mChildBuffer = populationBuffer + POPULATION_SIZE;
			mParentBuffer = parentBuffer;
			mSurvivorBuffer = survivorBuffer;

           // Initialise population
            mGenerationCount = 0;
            for(uint32_t i = 0; i < POPULATION_SIZE; ++i) {
				mPopulationBuffer[i] = population + i;
                Initialise(population[i].genome);
				population[i].isParent = 0;
				population[i].isChild = 0;
				population[i].fitness = CalculateFitness(population[i].genome);
				population[i].age = 0;
            }
			for(uint32_t i = 0; i < CHILDREN_PER_GENERATION; ++i) mChildBuffer[i] = children + i;

            while(! ShouldTerminate()) {
				OnGenerationBegin();

                // Generate children
				OnChildGeneration();

				// Select parents
				for(uint32_t i = 0; i < CHILDREN_PER_GENERATION; ++i) {
					SelectParents();
					for(int32_t j = 0; j < PARENTS_PER_CHILD; ++j) mParentBuffer[j]->isParent = 1;
					Crossover(children[i].genome);
					Mutate(children[i].genome);
					children[i].isParent = 0;
					children[i].isChild = 1;
					children[i].fitness = CalculateFitness(children[i].genome);
					children[i].age = 0;
				}

                // Select survivors
				OnSurvivorSelection();
                SelectSurvivors();
				population = population == mPopulationBufferA ? mPopulationBufferB : mPopulationBufferA;

                // End generation
				for(uint32_t i = 0; i < POPULATION_SIZE; ++i) {
					Phenotype* const phenotype = population + i;
					*phenotype = *survivorBuffer[i];
					populationBuffer[i] = phenotype;
					++phenotype->age;
					phenotype->isParent = 0;
					phenotype->isChild = 0;
				}
				OnGenerationEnd();
                ++mGenerationCount;
            }

			if(population == mPopulationBufferA) return mPopulationBufferA;
			return mPopulationBufferB;
        }
    };
}}

#endif
