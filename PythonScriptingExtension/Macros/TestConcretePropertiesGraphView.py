###################################
# This script creates all graph combinations for the concrete properties graph
###################################

import BridgeLink

pgsuper = BridgeLink.GetPGSuper()

graph_controller = pgsuper.Views.CreateGraphViewByName("Concrete Properties")
concrete_properties_graph = BridgeLink.ConcretePropertiesGraph.GetGraphController(graph_controller)

property_types = [BridgeLink.ConcreteProperty.Fc, BridgeLink.ConcreteProperty.Ec] #, BridgeLink.ConcreteProperty.Shrinkage, BridgeLink.ConcreteProperty.Creep (for time step analysis only)
element_types = [BridgeLink.ConcretePropertiesElement.Deck, BridgeLink.ConcretePropertiesElement.Segment] #, BridgeLink.ConcretePropertiesElement.Closure (for PGSplice only)
xaxis_types = [BridgeLink.ConcretePropertiesXAxis.TimeLinear,BridgeLink.ConcretePropertiesXAxis.TimeLog,BridgeLink.ConcretePropertiesXAxis.AgeLinear,BridgeLink.ConcretePropertiesXAxis.AgeLog,BridgeLink.ConcretePropertiesXAxis.Interval]

concrete_properties_graph.ShowGrid = False

nGroups = pgsuper.Bridge.GetGirderGroupCount();
for group in range(0,nGroups):
   nGirders = pgsuper.Bridge.GetGirderCount(group)
   for girder in range(0,nGirders):
      girderKey = BridgeLink.GirderKey(group,girder)
      nSegments = pgsuper.Bridge.GetSegmentCount(girderKey)
      for segment in range(0,nSegments):
         segmentKey = BridgeLink.SegmentKey(group,girder,segment)
         concrete_properties_graph.SelectSegment(segmentKey)
         for element in element_types:
            concrete_properties_graph.Element = element
            for property in property_types:
               concrete_properties_graph.Property = property
               for xaxis in xaxis_types:
                  concrete_properties_graph.XAxis = xaxis

concrete_properties_graph.ShowGrid = True
