#pragma once


///////////////////////////////////////////////////////////////////////////
// CIfcModelBuilder
class CIfcModelBuilder
{
public:
    CIfcModelBuilder(void);
    ~CIfcModelBuilder(void);

    enum SchemaType
    {
        Schema_4x3_rc3,
        //Schema_4x3_rc4
    };

    void BuildModel(IBroker* pBroker, const CString& strFilePath, SchemaType schemaType);

private:
    template <typename Schema>
    void BuildModel(IBroker* pBroker, const CString& strFilePath);
};

