###################################
# This script creates all graph combinations for the analysis results view
###################################


import BridgeLink

##################################################################
def CreateGraphs(girderKey):
   actionTypes = analysis_results_graph.GetActionTypes()
   analysis_types = [BridgeLink.AnalysisType.Simple, BridgeLink.AnalysisType.Continuous, BridgeLink.AnalysisType.Envelope]
   result_types = [BridgeLink.AnalysisResultsType.Cummulative, BridgeLink.AnalysisResultsType.Incremental]
   graphModes = [BridgeLink.AnalysisResultsGraphMode.Loading,BridgeLink.AnalysisResultsGraphMode.Interval]
   analysis_results_graph.SelectGirder(girderKey)
   for mode in graphModes: # iterate over each graph mode
      analysis_results_graph.GraphMode = mode
      for actionType in actionTypes: # iterate over each action type
         analysis_results_graph.ActionType = actionType
         graph_type_count = analysis_results_graph.GetGraphTypeCount()
         for graphType in range(0,graph_type_count): # iterate over each graph type for the current action type
            analysis_results_graph.SelectGraphTypeByIndex(graphType)
            graph_count = analysis_results_graph.GetGraphCount() # get the number of graphs for the current graph type
            graphs_to_select = list(range(0,graph_count)) # create a list of graphs to select (select all)
            analysis_results_graph.SelectGraphsByIndex(graphs_to_select) # select the graphs
            for analysis_type in analysis_types:
               analysis_results_graph.AnalysisType = analysis_type 
               for results_type in result_types:
                  analysis_results_graph.ResultsType = results_type
                  if actionType == BridgeLink.ActionType.Stress:
                     analysis_results_graph.SelectStressLocation(BridgeLink.StressLocation.TopGirder,False)
                     analysis_results_graph.SelectStressLocation(BridgeLink.StressLocation.BottomGirder,False)
                     analysis_results_graph.SelectStressLocation(BridgeLink.StressLocation.TopGirder,True)
                     analysis_results_graph.SelectStressLocation(BridgeLink.StressLocation.BottomGirder,True)
##################################################################

pgsuper = BridgeLink.GetPGSuper()

graph_controller = pgsuper.Views.CreateGraphViewByName("Analysis Results - On site")
analysis_results_graph = BridgeLink.AnalysisResultsGraph.GetGraphController(graph_controller)

## plot a specific graph
#girderKey = BridgeLink.GirderKey(BridgeLink.AllGroups(),0)
#analysis_results_graph.SelectGirder(girderKey)
#analysis_results_graph.GraphMode = BridgeLink.AnalysisResultsGraphMode.Interval
#analysis_results_graph.ActionType = BridgeLink.ActionType.Moment
#analysis_results_graph.SelectGraphTypeByName("Min Moment Capacity")
#graph_count = analysis_results_graph.GetGraphCount() # get the number of graphs for the current graph type
#graphs_to_select = list(range(0,graph_count)) # create a list of graphs to select (select all)
#analysis_results_graph.SelectGraphsByIndex(graphs_to_select) # select the graphs

# create all graphs for each girder line
nGirderLines = pgsuper.Bridge.GetGirderlineCount()
g = list()
g.append(0)
g.append(nGirderLines)
#for girder_line in range(0,nGirderLines):
for girder_line in g: # actually just doing first and last girder lines
   girderKey = BridgeLink.GirderKey(BridgeLink.AllGroups(),girder_line)
   CreateGraphs(girderKey)

# create all graphs for each individual girder
nGroups = pgsuper.Bridge.GetGirderGroupCount();
if 1 < nGroups: # if there is only one group, then no need to do this since it is exactly the same as above
    for group in range(0,nGroups):
       nGirders = pgsuper.Bridge.GetGirderCount(group)
       girders = list()
       girders.append(0)
       girders.append(nGirders-1)
       #for girder in range(0,nGirders):
       for girder in girders: # actually doing only first and last girder
             girderKey = BridgeLink.GirderKey(group,girder)
             CreateGraphs(girderKey)
