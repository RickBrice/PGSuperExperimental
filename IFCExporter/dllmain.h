// dllmain.h : Declaration of module class.

class CIFCExporterModule : public ATL::CAtlDllModuleT< CIFCExporterModule >
{
public :
	DECLARE_LIBID(LIBID_IFCExporterLib)
	DECLARE_REGISTRY_APPID_RESOURCEID(IDR_IFCEXPORTER, "{2084AE59-5080-4C70-B5A9-497113DEAEC6}")
};

extern class CIFCExporterModule _AtlModule;
