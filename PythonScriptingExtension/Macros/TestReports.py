import BridgeLink
pgsuper = BridgeLink.GetPGSuper()

nReports = pgsuper.Reports.Count()

for i in range(0,nReports):
   pgsuper.Views.CreateReport(i,False)



