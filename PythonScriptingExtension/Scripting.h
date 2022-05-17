#pragma once


// {83CABDBC-03B9-4991-AA5A-F8D67188B5DC}
DEFINE_GUID(IID_IScripting ,
   0x83cabdbc, 0x3b9, 0x4991, 0xaa, 0x5a, 0xf8, 0xd6, 0x71, 0x88, 0xb5, 0xdc);
interface IScripting : IUnknown
{
   virtual LPCTSTR GetScriptFolder() const = 0;
   virtual std::vector<CString> GetScriptNames() const = 0;
   virtual std::vector<CString> GetScriptFiles() const = 0;
   virtual CString GetScriptFile(LPCTSTR lpszScriptName) const = 0;
};