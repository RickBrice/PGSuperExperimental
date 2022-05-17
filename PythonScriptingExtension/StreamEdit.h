#pragma once
#include <iosfwd>
#include <streambuf>

// CStreamEdit

class CStreamEdit : public CEdit, public std::streambuf
{
	DECLARE_DYNAMIC(CStreamEdit)

public:
	CStreamEdit();
	virtual ~CStreamEdit();

   virtual std::streamsize xsputn(const char *_Ptr, std::streamsize _Count) override;

protected:
	DECLARE_MESSAGE_MAP()
};


