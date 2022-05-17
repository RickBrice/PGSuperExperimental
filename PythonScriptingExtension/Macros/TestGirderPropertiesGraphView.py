###################################
# This script creates all graph combinations for the girder properties graph
###################################

import BridgeLink

pgsuper = BridgeLink.GetPGSuper()

graph_controller = pgsuper.Views.CreateGraphViewByName("Girder Properties")
girder_properties_graph = BridgeLink.GirderPropertiesGraph.GetGraphController(graph_controller)

girder_properties_graph.ShowGrid = False
girder_properties_graph.ShowGirder = False

firstInterval = girder_properties_graph.GetFirstInterval()
lastInterval = girder_properties_graph.GetLastInterval()
intervals = list(range(firstInterval,lastInterval+1))

property_types = [BridgeLink.GirderPropertyType.Height, \
    BridgeLink.GirderPropertyType.Area,\
    BridgeLink.GirderPropertyType.MomentOfInertia,\
    BridgeLink.GirderPropertyType.Centroid,\
    BridgeLink.GirderPropertyType.SectionModulus,\
    BridgeLink.GirderPropertyType.KernPoint,\
    BridgeLink.GirderPropertyType.StrandEccentricity,\
    #BridgeLink.GirderPropertyType.TendonEccentricity,\
    #BridgeLink.GirderPropertyType.TendonProfile,\
    BridgeLink.GirderPropertyType.EffectiveFlangeWidth,\
    BridgeLink.GirderPropertyType.Fc,\
    BridgeLink.GirderPropertyType.Ec]

section_property_types = [BridgeLink.SectionPropertyType.Gross, \
    BridgeLink.SectionPropertyType.Transformed, \
    BridgeLink.SectionPropertyType.NetGirder, \
    BridgeLink.SectionPropertyType.NetDeck]

nGroups = pgsuper.Bridge.GetGirderGroupCount();
for group in range(0,nGroups):
   nGirders = pgsuper.Bridge.GetGirderCount(group)
   for girder in range(0,nGirders):
         girderKey = BridgeLink.GirderKey(group,girder)
         girder_properties_graph.SelectGirder(girderKey)
         for interval in intervals:
            girder_properties_graph.Interval = interval
            for property in property_types:
               girder_properties_graph.GirderProperty = property
               if girder_properties_graph.IsInvariantProperty(property) != True:
                  for section_property in section_property_types:
                     girder_properties_graph.SectionProperty = section_property


girder_properties_graph.ShowGrid = True
girder_properties_graph.ShowGirder = True
