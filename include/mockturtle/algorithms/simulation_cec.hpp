/* mockturtle: C++ logic network library
* Copyright (C) 2018-2021  EPFL
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
* OTHER DEALINGS IN THE SOFTWARE.
*/

/*!
 \file simulation_cec.hpp
 \brief Simulation-based CEC

 EPFL CS-472 2021 Final Project Option 2
*/

#pragma once

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../utils/node_map.hpp"
#include "miter.hpp"
#include "simulation.hpp"
namespace mockturtle
{

/* Statistics to be reported */
    struct simulation_cec_stats
    {
        /*! \brief Split variable (simulation size). */
        uint32_t split_var{ 0 };

        /*! \brief Number of simulation round. */
        uint32_t rounds{ 0 };
    };

    namespace detail
    {

        class OptimizedSimulator {
        public:
            unsigned int numVars;
            unsigned int splitVar;
            unsigned int round;

            OptimizedSimulator( unsigned int numVars, unsigned int splitVar, unsigned int round )
                    : numVars( numVars ), splitVar(splitVar), round( round ) {}

            kitty::dynamic_truth_table compute_constant( bool value ) const
            {
                kitty::dynamic_truth_table tt( numVars );
                return value ? ~tt : tt;
            }

            kitty::dynamic_truth_table compute_pi( uint32_t index ) const
            {
                kitty::dynamic_truth_table tt(numVars);
                if(index < splitVar) {
                    kitty::create_nth_var(tt, index);
                } else {
                    tt = getBitAtIndex(index) ? tt : ~tt;
                }
                return tt;
            }

            kitty::dynamic_truth_table compute_not( kitty::dynamic_truth_table const& value ) const
            {
                return ~value;
            }

        private:
            // returns true
            bool getBitAtIndex(unsigned int index) const {
                int shifted = round >> (index - splitVar);
                int bit = shifted & 1;
                return bit;
            }
        };

        template<class Ntk>
        class simulation_cec_impl
        {
        public:
            using pattern_t = unordered_node_map<kitty::dynamic_truth_table, Ntk>;
            using node = typename Ntk::node;
            using signal = typename Ntk::signal;

        public:
            explicit simulation_cec_impl( Ntk& ntk, simulation_cec_stats& st )
                    : _ntk( ntk ),
                      _st( st )
            {
            }

            bool run()
            {
                bool result = false;
                // Initial parameter calculation
                calculate_st();
                // Pattern generation - Cofactoring - then simulating the smaller truth tables
                result = sim();
                return result;
            }

        private:
            void calculate_st() {
                // Variables stored in _st
                // Calculating the splitting variable and required number of round to avoid a large simulation
                uint64_t  n = _ntk.num_pis();
                if (n <= 6)
                    _st.split_var = n;
                else {
                    uint64_t m,V;
                    V = _ntk.size();
                    m = 7;
                    while (m < n && (32 + pow(2,(m-3 + 1 ))) * V <= pow(2, 29)) {   //Fixed bug, you need to anticipate by 1 cycle or else it gets m +1 instead of m
                        m++;
                    }
                    _st.split_var = m;
                }

                _st.rounds = pow(2,n - _st.split_var);
            }
            bool sim() {
                // Define the truthtable that we have
                for(int round = 0; round < _st.rounds; round++){
                    OptimizedSimulator optimizedSimulator(_ntk.num_pis(), _st.split_var, round);
                    const std::vector<kitty::dynamic_truth_table> tts = simulate<kitty::dynamic_truth_table>(_ntk, optimizedSimulator);
                    for(auto& po : tts) {
                        if(!kitty::is_const0(po)){
                            return false;
                        }
                    }
                }
                return true;
            }

        private:
            Ntk& _ntk;
            simulation_cec_stats& _st;
            /* you can add other attributes here */
        };

    } // namespace detail

/* Entry point for users to call */

/*! \brief Simulation-based CEC.
*
* This function implements a simulation-based combinational equivalence checker.
* The implementation creates a miter network and run several round of simulation
* to verify the functional equivalence. For memory and speed reasons this approach
* is limited up to 40 input networks. It returns an optional which is `nullopt`,
* if the network has more than 40 inputs.
*/
    template<class Ntk>
    std::optional<bool> simulation_cec( Ntk const& ntk1, Ntk const& ntk2, simulation_cec_stats* pst = nullptr )
    {
        static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
        static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
        static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
        static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
        static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
        static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
        static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
        static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

        simulation_cec_stats st;

        bool result = false;

        if ( ntk1.num_pis() > 40 )
            return std::nullopt;

        auto ntk_miter = miter<Ntk>( ntk1, ntk2 );

        if ( ntk_miter.has_value() )
        {
            detail::simulation_cec_impl p( *ntk_miter, st );
            result = p.run();
        }

        if ( pst )
            *pst = st;

        return result;
    }

} // namespace mockturtle
