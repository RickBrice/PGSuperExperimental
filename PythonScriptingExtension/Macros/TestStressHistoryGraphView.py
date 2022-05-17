###################################
# This script creates all graph combinations for the stress history graph
###################################

import BridgeLink

pgsuper = BridgeLink.GetPGSuper()

graph_controller = pgsuper.Views.CreateGraphViewByName("Stress History")
stress_history_graph = BridgeLink.StressHistoryGraph.GetGraphController(graph_controller)
xaxis_types = [BridgeLink.StressHistoryXAxis.TimeLinear,BridgeLink.StressHistoryXAxis.TimeLog,BridgeLink.StressHistoryXAxis.Interval]

stress_history_graph.ShowGrid = False

nSpans = pgsuper.Bridge.GetSpanCount();
for span in range(0,nSpans):
   nGirders = pgsuper.Bridge.GetGirderCountBySpan(span)
   for girder in range(0,nGirders):
      spanKey = BridgeLink.SpanKey(span,girder)
      points_of_interest = pgsuper.PointOfInterest.GetSpanPointsOfInterest(spanKey)
      for poi in points_of_interest:
         stress_history_graph.SelectLocation(poi)
         for xaxis in xaxis_types:
           stress_history_graph.XAxis = xaxis
           stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.TopDeck,True)
           stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.BottomDeck,True)
           stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.TopGirder,True)
           stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.BottomGirder,True)
           stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.TopDeck,False)
           stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.BottomDeck,False)
           stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.TopGirder,False)
           stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.BottomGirder,False)

stress_history_graph.ShowGrid = True
