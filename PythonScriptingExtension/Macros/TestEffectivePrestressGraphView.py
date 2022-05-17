###################################
# This script creates all graph combinations for the effective prestress graph
###################################

import BridgeLink

pgsuper = BridgeLink.GetPGSuper()

graph_controller = pgsuper.Views.CreateGraphViewByName("Effective Prestress")
effective_prestress_graph = BridgeLink.EffectivePrestressGraph.GetGraphController(graph_controller)

graph_modes = [BridgeLink.EffectivePrestressGraphMode.Force,BridgeLink.EffectivePrestressGraphMode.Stress]
strand_types = [BridgeLink.EffectivePrestressStrandType.Permanent,BridgeLink.EffectivePrestressStrandType.Temporary]

minInterval = effective_prestress_graph.GetMinimumInterval()
maxInterval = effective_prestress_graph.GetMaximumInterval()
intervals = list(range(minInterval,maxInterval+1))

effective_prestress_graph.SelectIntervals(intervals)

effective_prestress_graph.ShowGrid = False
effective_prestress_graph.ShowGirder = False

nGroups = pgsuper.Bridge.GetGirderGroupCount();
for group in range(0,nGroups):
   nGirders = pgsuper.Bridge.GetGirderCount(group)
   for girder in range(0,nGirders):
         girderKey = BridgeLink.GirderKey(group,girder)
         effective_prestress_graph.SelectGirder(girderKey)
         for strand in strand_types:
            effective_prestress_graph.StrandType = strand
            for mode in graph_modes:
               effective_prestress_graph.Mode = mode

effective_prestress_graph.ShowGrid = True
effective_prestress_graph.ShowGirder = True
