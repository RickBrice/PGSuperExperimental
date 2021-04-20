// dllmain.cpp : Implementation of DllMain.

#include "stdafx.h"
#include "resource.h"
#include "IFCExporter_i.h"
#include "dllmain.h"

CIFCExporterModule _AtlModule;

class CIFCExporterApp : public CWinApp
{
public:

// Overrides
	virtual BOOL InitInstance();
	virtual int ExitInstance();

	DECLARE_MESSAGE_MAP()
};

BEGIN_MESSAGE_MAP(CIFCExporterApp, CWinApp)
END_MESSAGE_MAP()

CIFCExporterApp theApp;

BOOL CIFCExporterApp::InitInstance()
{
	return CWinApp::InitInstance();
}

int CIFCExporterApp::ExitInstance()
{
	return CWinApp::ExitInstance();
}
