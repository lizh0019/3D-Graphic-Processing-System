// openglpracticeDoc.cpp : implementation of the COpenglpracticeDoc class
//

#include "stdafx.h"
#include "openglpractice.h"

#include "openglpracticeDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeDoc

IMPLEMENT_DYNCREATE(COpenglpracticeDoc, CDocument)

BEGIN_MESSAGE_MAP(COpenglpracticeDoc, CDocument)
	//{{AFX_MSG_MAP(COpenglpracticeDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeDoc construction/destruction

COpenglpracticeDoc::COpenglpracticeDoc()
{
	// TODO: add one-time construction code here

}

COpenglpracticeDoc::~COpenglpracticeDoc()
{
}

BOOL COpenglpracticeDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}



/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeDoc serialization

void COpenglpracticeDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeDoc diagnostics

#ifdef _DEBUG
void COpenglpracticeDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void COpenglpracticeDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeDoc commands
