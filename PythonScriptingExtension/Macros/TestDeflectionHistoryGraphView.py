###################################
# This script creates all graph combinations for the deflection history graph
###################################

import BridgeLink

pgsuper = BridgeLink.GetPGSuper()

graph_controller = pgsuper.Views.CreateGraphViewByName("Deflection History")
deflection_history_graph = BridgeLink.DeflectionHistoryGraph.GetGraphController(graph_controller)
xaxis_types = [BridgeLink.DeflectionHistoryXAxis.TimeLinear,BridgeLink.DeflectionHistoryXAxis.TimeLog,BridgeLink.DeflectionHistoryXAxis.Interval]

deflection_history_graph.ShowGrid = False

nSpans = pgsuper.Bridge.GetSpanCount();
for span in range(0,nSpans):
   nGirders = pgsuper.Bridge.GetGirderCountBySpan(span)
   for girder in range(0,nGirders):
      spanKey = BridgeLink.SpanKey(span,girder)
      points_of_interest = pgsuper.PointOfInterest.GetSpanPointsOfInterest(spanKey)
      for poi in points_of_interest:
         deflection_history_graph.SelectLocation(poi)
         for xaxis in xaxis_types:
           deflection_history_graph.XAxis = xaxis

deflection_history_graph.ShowGrid = True
