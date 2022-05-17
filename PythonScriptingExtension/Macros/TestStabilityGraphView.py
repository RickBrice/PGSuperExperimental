###################################
# This script creates all graph combinations for the concrete properties graph
###################################

import BridgeLink

pgsuper = BridgeLink.GetPGSuper()

graph_controller = pgsuper.Views.CreateGraphViewByName("Girder Stability")
stability_graph = BridgeLink.StabilityGraph.GetGraphController(graph_controller)

graph_modes = [BridgeLink.StabilityGraphMode.Lifting, BridgeLink.StabilityGraphMode.Hauling]

stability_graph.ShowGrid = False

nGroups = pgsuper.Bridge.GetGirderGroupCount();
for group in range(0,nGroups):
   nGirders = pgsuper.Bridge.GetGirderCount(group)
   for girder in range(0,nGirders):
      girderKey = BridgeLink.GirderKey(group,girder)
      nSegments = pgsuper.Bridge.GetSegmentCount(girderKey)
      for segment in range(0,nSegments):
         segmentKey = BridgeLink.SegmentKey(group,girder,segment)
         stability_graph.SelectSegment(segmentKey)
         for mode in graph_modes:
            stability_graph.Mode = mode

stability_graph.ShowGrid = True
