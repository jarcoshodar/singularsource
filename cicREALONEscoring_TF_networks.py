#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 11:25:19 2023

@author: jarcos
"""

import os
import pandas as pd
import networkx as nx
import re


rootDir = '/home/CICBIOGUNE/jarcos/Documentos/Real_Data_2/Rapamycin/results3/Integrated'

results_root_directory = os.path.expanduser('~/TFscoresdel/')
masteregulator_directory = os.path.expanduser('~/MRflipsdel')
os.makedirs(results_root_directory, exist_ok=True)
os.makedirs(masteregulator_directory, exist_ok=True)

# Walk through the directory
for dirName, subdirList, fileList in os.walk(rootDir):
    if not subdirList and 'GRNI5' in dirName:
        parts = dirName.split('/')
        organ = parts[-2]
        study = parts[-5]
        GRNI = parts[-1]

        # Construct the results directory by appending the subfolder, GRNI and tissue type
        results_directory = os.path.join(results_root_directory, study, GRNI)

        os.makedirs(results_directory, exist_ok=True)

        # Loop over all CSV files in the directory
        for filename in fileList:
            if not filename.endswith('.csv') or 'FAIL' in filename:
                continue  # Skip non-CSV files
            filepath = os.path.join(dirName, filename)
            df = pd.read_csv(filepath)
            if len(df) == 0:
               continue
            #we initialize the network as object G, a directed graph
            G = nx.DiGraph()
            celltype = re.search('GRNI5([A-Z][a-z]*)', filename).group(1)
            #unlike in prune_gurobi we value them equally and just try
            #to keep track of active inhibitions and give them priority later
            effect_mapping = {'Activation': 1, 'Inhibition': -1}
            
            #fill the network with the dataframe from GRNOpt
            #based on column names
            for _, row in df.iterrows():
                G.add_edge(row['TF'], row['target'], weight=effect_mapping[row['effect']])
            
            #my default they are all neutral, or name also for undefined.
            #we activate manually or let the network dictate activations
            node_states = {node: 'neutral' for node in G.nodes()}
            # store TF/regulon scores in a ditcionary
            node_scores = {}
            
            
            def determine_state(node):
                """
                Determines the next state of a given node based on the current state of its parent nodes.
                
                The function checks all edges coming into the node. If any of these edges come from an 'active'
                node and have a negative weight (representing inhibition), the node will become 'inactive'.
                This priority is given to replicate the 'Inhibition Dominance' rule of GRNOpt.
                If any edges come from an 'active' node and have a positive weight (representing activation),
                and there are no active inhibitions, the node will become 'active'. If there are no active 
                influences at all, the node retains its current state, neutral by design of the script
            
                Parameters:
                node (str): The node for which the next state is to be determined.
            
                Returns:
                str: The next state of the node ('active', 'inactive', or its current state).
                """
                parent_edges = G.in_edges(node, data=True)
                active_activations = [edge_data['weight'] > 0 for u, _, edge_data in parent_edges if node_states[u] == 'active']
                active_inhibitions = [edge_data['weight'] < 0 for u, _, edge_data in parent_edges if node_states[u] == 'active']
            
                if any(active_inhibitions):
                    return 'inactive'
                elif any(active_activations):
                    return 'active'
                else:
                    return node_states[node]  # Returns current state if not certain
            
            def hash_node_states(node_states):
                """
                Creates a hashable representation of the current states of all nodes in the network.
                
                The function uses a frozenset, which is an immutable and hashable version of a Python set.
                This allows us to use the set of current states as a dictionary key or in a set.
            
                Parameters:
                node_states (dict): A dictionary mapping each node to its current state.
            
                Returns:
                int: A hashable representation of the current states of all nodes.
                """
                return hash(frozenset(node_states.items()))
            
            def dfs_update_states(start_node, visited=None):
                """
                Performs a depth-first search (DFS) through the network to update the states of all nodes,
                starting from a given node.
                """
                global node_states
                if visited is None:
                    visited = set()
                visited.add(start_node)
                
                for _, child_node, edge_data in G.out_edges(start_node, data=True):
                    if child_node not in visited:
                        node_states[child_node] = determine_state(child_node)
                        dfs_update_states(child_node, visited)
        
        
            # Initialize dict to store states. some genes require several steps
            #to do all the updates, so recall to always refer to a network as
            #all_states['gene name'][-1]
            #so that a unique 1 index and the correct 3 index are correctly called
            all_states = {}
            
            # Load the boolean data here
            def grab_bool(old_filename):
                old_dir, old_file = os.path.split(old_filename)
        
                if old_file.startswith('GRNI5'):
                    new_file = old_file[5:]
                else:
                    new_file = old_file
        
                new_file = 'bool_' + new_file
        
                new_dir = os.path.dirname(old_dir)
        
                new_dir = os.path.join(new_dir, 'booldata')
        
                new_filename = os.path.join(new_dir, new_file)
        
                return new_filename
        
            bool_grni = grab_bool(filepath)
            booldata = pd.read_csv(bool_grni)
        
            # Create a dictionary to map gene names to their boolean states
            # Only include genes that are also nodes in the network graph
            bool_states = {row[0]: row[1] for _, row in booldata.iterrows() if row[0] in G.nodes()}
        
            all_states = {}
        
            #We iterate for every node in the graph
            for node in G.nodes:
                #set the node to active and track it's consequences
                node_states[node] = 'active'
        
                #set of processed outcomes
                seen_states = set()
        
                #initialize list to store states for this node
                node_states_list = []
        
                #repeat until the states of the nodes no longer change
                previous_states = {}
                #this must be done to avoid infinite recursion due to cyclical gene networks
                while hash_node_states(previous_states) != hash_node_states(node_states):
                    #check if current state has been seen before
                    current_state_hash = hash_node_states(node_states)
                    if current_state_hash in seen_states:
                        print(f"Cycle detected after activating node {node}. Breaking loop.")
                        break
        
                    #store the state in the while
                    node_states_list.append(node_states.copy())
        
                    seen_states.add(current_state_hash) #hash for future reference
                    previous_states = node_states.copy() #.copy() avoids modification on future redefinition across loops
                    dfs_update_states(node) #implementation of depth first search goes here
        
                #state list for this node. as mentioned before, always refer to last for real results
                all_states[node] = node_states_list
        
                # Here's where we update the scoring
                # Score is the number of nodes that match the boolean data
                # print(f"Processing node: {node}")
                # print(f"Node states: {node_states}")
                # print(f"Boolean states: {bool_states}")
                # scored_genes = set()
                # score = 0
                # for gene, state in node_states.items():
                #     if gene in bool_states and ((state == 'active' and bool_states[gene] == 1) or (state == 'inactive' and bool_states[gene] == 0)):
                #         if gene not in scored_genes:
                #             score += 1
                #             scored_genes.add(gene)
                # score /= len(node_states)
        
                score = sum(1 for gene, state in node_states.items() if gene in bool_states and ((state == 'active' and bool_states[gene] == 1) or (state == 'inactive' and bool_states[gene] == 0))) / len(node_states)
                # print(f"Score for node {node}: {score}")
        
                node_scores[node] = score
                # score = sum(1 for gene, state in node_states.items() if gene in bool_states and ((state == 'active' and bool_states[gene] == 1) or (state != 'active' and bool_states[gene] == 0)))
                # node_scores[node] = score
        
                #reset all nodes back to 'neutral' before moving on to the next node
                node_states = {node: 'neutral' for node in G.nodes()}
            
            # Now let's compare
            booldata = pd.read_csv(bool_grni)
            # Convert the node_scores dict to a dataframe
            node_scores_df = pd.DataFrame.from_dict(node_scores, orient='index', columns=['Score'])
        
            # Sort the DataFrame by score
            node_scores_df = node_scores_df.sort_values(by='Score', ascending=False)
            
            # Add the 'Orig' column
            node_scores_df['Orig'] = [bool_states.get(tf, 'N/A') for tf in node_scores_df.index]
        
            # Save the DataFrame to a CSV file
            output_filepath = os.path.join(results_directory, f'{filename}_scores.csv')
            node_scores_df.to_csv(output_filepath, sep='\t', index=True, index_label='TF')
            
            
            df2 = pd.DataFrame(list(node_scores.items()), columns=['TF', 'Score'])

            #Merge dataframes on 'TF'
            merged = pd.merge(df, df2, on='TF', how='left')
            print(merged)
            print(filepath)
            # Save the merged dataframe
            merged.to_csv(filepath, index=False)
                        
            max_score = node_scores_df['Score'].max()
        
            # Get all top-scoring TFs (there could be a tie)
            top_scoring_TFs = node_scores_df[node_scores_df['Score'] == max_score].index.tolist()
        
            # For each top-scoring TF, activate it and save the genes that are not 'neutral'
            for tf in top_scoring_TFs:
                # Activate the TF
                node_states[tf] = 'active'
                # Run the DFS to update states
                dfs_update_states(tf)
                # Find all non-neutral genes
                non_neutral_genes = [node for node, state in node_states.items() if state != 'neutral']
                # Reset the states for the next TF
                node_states = {node: 'neutral' for node in G.nodes()}
                # Save the non-neutral genes to a file
                output_filename = f"{study}_{organ}_{celltype}_{tf}_non_neutral_genes.txt"

                output_filepath = os.path.join(masteregulator_directory, output_filename)
                
                with open(output_filepath, "w") as f:
                    for gene in non_neutral_genes:
                        f.write(f"{gene}\n")
