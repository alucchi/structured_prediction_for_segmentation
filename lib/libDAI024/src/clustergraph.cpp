/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <set>
#include <vector>
#include <iostream>
#include <dai/varset.h>
#include <dai/clustergraph.h>


namespace dai {


using namespace std;


ClusterGraph::ClusterGraph( const std::vector<VarSet> & cls ) : G(), vars(), clusters() {
    // construct vars, clusters and edge list
    vector<Edge> edges;
    daiforeach( const VarSet &cl, cls ) {
        if( find( clusters.begin(), clusters.end(), cl ) == clusters.end() ) {
            // add cluster
            size_t n2 = clusters.size();
            clusters.push_back( cl );
            for( VarSet::const_iterator n = cl.begin(); n != cl.end(); n++ ) {
                size_t n1 = find( vars.begin(), vars.end(), *n ) - vars.begin();
                if( n1 == vars.size() )
                    // add variable
                    vars.push_back( *n );
                edges.push_back( Edge( n1, n2 ) );
            }
        } // disregard duplicate clusters
    }

    // Create bipartite graph
    G.construct( vars.size(), clusters.size(), edges.begin(), edges.end() );
}


size_t ClusterGraph::eliminationCost( size_t i ) const {
    return eliminationCost_MinFill( *this, i );
}


ClusterGraph ClusterGraph::VarElim( const std::vector<Var> &ElimSeq ) const {
    return VarElim( sequentialVariableElimination( ElimSeq ) );
}


ClusterGraph ClusterGraph::VarElim_MinFill() const {
    return VarElim( greedyVariableElimination( &eliminationCost_MinFill ) );
}


size_t sequentialVariableElimination::operator()( const ClusterGraph &cl, const std::set<size_t> &/*remainingVars*/ ) {
    return cl.findVar( seq.at(i++) );
}


size_t greedyVariableElimination::operator()( const ClusterGraph &cl, const std::set<size_t> &remainingVars ) {
    set<size_t>::const_iterator lowest = remainingVars.end();
    size_t lowest_cost = -1UL;
    for( set<size_t>::const_iterator i = remainingVars.begin(); i != remainingVars.end(); i++ ) {
        size_t cost = heuristic( cl, *i );
        if( lowest == remainingVars.end() || lowest_cost > cost ) {
            lowest = i;
            lowest_cost = cost;
        }
    }
    return *lowest;
}


size_t eliminationCost_MinNeighbors( const ClusterGraph &cl, size_t i ) {
    std::vector<size_t> id_n = cl.G.delta1( i );
    return id_n.size();
}


size_t eliminationCost_MinWeight( const ClusterGraph &cl, size_t i ) {
    std::vector<size_t> id_n = cl.G.delta1( i );
    
    size_t cost = 1;
    for( size_t _i = 0; _i < id_n.size(); _i++ )
        cost *= cl.vars[id_n[_i]].states();

    return cost;
}


size_t eliminationCost_MinFill( const ClusterGraph &cl, size_t i ) {
    std::vector<size_t> id_n = cl.G.delta1( i );

    size_t cost = 0;
    // for each unordered pair {i1,i2} adjacent to n
    for( size_t _i1 = 0; _i1 < id_n.size(); _i1++ )
        for( size_t _i2 = _i1 + 1; _i2 < id_n.size(); _i2++ ) {
            // if i1 and i2 are not adjacent, eliminating n would make them adjacent
            if( !cl.adj(id_n[_i1], id_n[_i2]) )
                cost++;
        }

    return cost;
}


size_t eliminationCost_WeightedMinFill( const ClusterGraph &cl, size_t i ) {
    std::vector<size_t> id_n = cl.G.delta1( i );

    size_t cost = 0;
    // for each unordered pair {i1,i2} adjacent to n
    for( size_t _i1 = 0; _i1 < id_n.size(); _i1++ )
        for( size_t _i2 = _i1 + 1; _i2 < id_n.size(); _i2++ ) {
            // if i1 and i2 are not adjacent, eliminating n would make them adjacent
            if( !cl.adj(id_n[_i1], id_n[_i2]) )
                cost += cl.vars[id_n[_i1]].states() * cl.vars[id_n[_i2]].states();
        }

    return cost;
}


} // end of namespace dai
