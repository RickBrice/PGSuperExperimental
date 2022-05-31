# MoveBridgeModelSectionCut
# This method moves the bridge model view section cut through a series
# of locations, starting at start, ending at end, and moving inc each time.
# This method was developed to show an example of using user defined Python methods
# in a macro
def MoveBridgeModelSectionCut(bridge_model_view,start,end,inc) :
   for i in range(start,end,inc) :
      bridge_model_view.CutStation = i


# This is where the macro starts running
if __name__ == "__main__" :
   import BridgeLink # import the BridgeLink module

   # get the pgsuper object
   pgsuper = BridgeLink.GetPGSuper()

   print(pgsuper.VersionString)
   print(pgsuper.VersionStringWithBuildNumber)
   print(pgsuper.Version)
   print(pgsuper.VersionWithBuildNumber)

   # get the bridge model view and move the section cut by
   # calling a user defined Python method
   bridge_model_view = pgsuper.Views.CreateBridgeModelView()
   MoveBridgeModelSectionCut(bridge_model_view,100,200,10)

   bridge_model_view.SelectDeck()

   # Edit the project properties
   pgsuper.Properties.BridgeName = "MyBridge"

   # Edit the environment settings
   pgsuper.Environment.ExposureCondition = BridgeLink.ExposureCondition.Severe
   pgsuper.Environment.RelHumidity = 85

   # set load modifiers
   pgsuper.LoadModifiers.SetDuctilityModifier(BridgeLink.LoadModifierLevel.Low,0.95)
   pgsuper.LoadModifiers.SetImportanceModifier(BridgeLink.LoadModifierLevel.Normal,0.99)
   pgsuper.LoadModifiers.SetRedundancyModifier(BridgeLink.LoadModifierLevel.High,1.05)

   print("Ductity load modifiers")
   print(pgsuper.LoadModifiers.GetDuctilityLevel())
   print(pgsuper.LoadModifiers.GetDuctilityFactor())

   # define the alignment
   alignment = BridgeLink.Alignment() # create a new alignment object
   alignment.RefStation = 500 # set the reference station

   # create an define a horizontal curve
   hc = BridgeLink.HorzCurve()
   hc.PIStation = 100
   hc.Radius = 500
   hc.FwdTangent = 0.3253

   # add the horizontal curve to the alignment
   alignment.HorzCurves.append(hc)

   # define the profile
   profile = BridgeLink.Profile()
   profile.Station = 100
   profile.Elevation = 40
   profile.Grade = 2.345/100.
   vc = BridgeLink.VertCurve()
   vc.PVIStation = 100
   vc.ExitGrade = -3.453/100.
   vc.L1 = 400
   profile.VertCurves.append(vc)

   # define superelevations
   roadwaySections = BridgeLink.RoadwaySections()
   roadwaySections.SegmentsPerSection = 2
   roadwaySections.ControllingRidgePoint = 1
   section = BridgeLink.SectionTemplate()
   section.Station = 100
   section.LeftSlope = -0.02
   section.RightSlope = -0.02
   roadwaySections.Templates.append(section)
   
   section.Station = 200
   section.LeftSlope = -0.02
   section.RightSlope = 0.02
   roadwaySections.Templates.append(section)

   section.Station = 300
   section.LeftSlope = -0.06
   section.RightSlope = 0.06
   roadwaySections.Templates.append(section)

   # assign the roadway elements to the project
   pgsuper.Roadway.Alignment = alignment
   pgsuper.Roadway.Profile = profile
   pgsuper.Roadway.RoadwaySections = roadwaySections

   # change the view mode of the bridge model view
   bridge_model_view.ViewMode = BridgeLink.BridgeModelViewMode.Alignment

   # change the bridge model view settings
   bridge_model_view.ShowBridge = True
   bridge_model_view.Schematic = False

   girder_model_view = pgsuper.Views.CreateGirderView()
   girder_model_view.Event = 3
   girder_model_view.Schematic = True
   girder_model_view.ShowCG = False
   girder_model_view.CutAtNext()
   girder_model_view.CutAtNext()
   girder_model_view.CutAtNext()

   loads_view = pgsuper.Views.CreateLoadsView()

   graph_controller = pgsuper.Views.CreateGraphViewByName("Girder Stability")
   stability_graph = BridgeLink.StabilityGraph.GetGraphController(graph_controller)
   stability_graph.Mode = BridgeLink.StabilityGraphMode.Hauling
   stability_graph.ShowGrid = False
   
   # stability_graph and graph_controller are the same object
   stability_graph.Minimize
   graph_controller.Maximize

   # effective prestress graph
   graph_controller = pgsuper.Views.CreateGraphViewByName("Effective Prestress")
   effective_prestress_graph = BridgeLink.EffectivePrestressGraph.GetGraphController(graph_controller)
   effective_prestress_graph.Mode = BridgeLink.EffectivePrestressGraphMode.Force

   minInterval = effective_prestress_graph.GetMinimumInterval()
   maxInterval = effective_prestress_graph.GetMaximumInterval()
   intervals = list(range(minInterval,maxInterval+1))

   effective_prestress_graph.SelectIntervals(intervals)

   graph_controller = pgsuper.Views.CreateGraphViewByName("Concrete Properties")
   concrete_properties_graph = BridgeLink.ConcretePropertiesGraph.GetGraphController(graph_controller)
   concrete_properties_graph.XAxis = BridgeLink.ConcretePropertiesXAxis.TimeLinear
   
   segmentKey = BridgeLink.SegmentKey(0,1,0)
   concrete_properties_graph.SelectSegment(segmentKey)
   concrete_properties_graph.Property = BridgeLink.ConcreteProperty.Ec
   concrete_properties_graph.ShowGrid = False

   # write test code for deflection and stress history after we get the poi thing worked out
   # graph_controller = pgsuper.Views.CreateGraphViewByName("Stress History")
   # stress_history_graph = BridgeLink.StressHistoryGraph.GetGraphController(graph_controller)
   # stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.TopGirder,True)
   # stress_history_graph.SelectStressLocation(BridgeLink.StressLocation.BottomGirder,True)

   graph_controller = pgsuper.Views.CreateGraphViewByName("Analysis Results - At Casting Yard")
   analysis_results_graph = BridgeLink.SegmentAnalysisResultsGraph.GetGraphController(graph_controller)
   analysis_results_graph.ShowGrid = False
   analysis_results_graph.ShowGirder = False
   analysis_results_graph.ResultsType = BridgeLink.AnalysisResultsType.Incremental
   analysis_results_graph.AnalysisType = BridgeLink.AnalysisType.Simple
   analysis_results_graph.SelectStressLocation(BridgeLink.StressLocation.TopGirder,True)
   analysis_results_graph.SelectStressLocation(BridgeLink.StressLocation.BottomGirder,False)
   analysis_results_graph.ActionType = BridgeLink.ActionType.Stress
   
   analysis_results_graph.GraphMode = BridgeLink.AnalysisResultsGraphMode.Interval # plot by interval
   graph_count = analysis_results_graph.GetGraphCount() # get the number of graphs
   graphs_to_select = list(range(0,graph_count)) # create a list of graphs to select
   analysis_results_graph.SelectGraphsByIndex(graphs_to_select) # select the graphs

   graph_names = list()
   selected_graphs = analysis_results_graph.GetSelectedGraphs()
   for i in selected_graphs:
      graph_names.append(analysis_results_graph.GetGraphName(i))
      msg = "Graph " + str(i) + " " + graph_names[i]
      print(msg)

   print(graph_names)

   analysis_results_graph.SelectGraphByIndex(0) # this doesn't seem to work

   # select every other graph by name
   graph_names = list()
   for i in range(0,graph_count,2):
      graph_names.append(analysis_results_graph.GetGraphName(i))

   analysis_results_graph.SelectGraphsByName(graph_names)
