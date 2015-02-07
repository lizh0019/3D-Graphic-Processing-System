// kettleDoc.cpp : implementation of the CKettleDoc class
//

#include "stdafx.h"
#include "kettle.h"
// Standard include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#include "kettleDoc.h"
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CKettleDoc

IMPLEMENT_DYNCREATE(CKettleDoc, CDocument)

BEGIN_MESSAGE_MAP(CKettleDoc, CDocument)
	//{{AFX_MSG_MAP(CKettleDoc)
	ON_COMMAND(ID_FILE_SAVE, OnSave3drma)
	ON_COMMAND(ID_OPEN_3DRMA, OnOpen3drma)
	ON_COMMAND(ID_FILE_OPEN, OnOpen3drma)
	ON_COMMAND(IDM_PROCESS, OnProcess)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CKettleDoc construction/destruction

CKettleDoc::CKettleDoc()
{
	// TODO: add one-time construction code here
opensuccess=FALSE;
}

CKettleDoc::~CKettleDoc()
{
}

BOOL CKettleDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)


	return TRUE;
}



/////////////////////////////////////////////////////////////////////////////
// CKettleDoc serialization

void CKettleDoc::Serialize(CArchive& ar)
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
// CKettleDoc diagnostics

#ifdef _DEBUG
void CKettleDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CKettleDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CKettleDoc commands


	/************************************************************************/
	/* Æ½ÒÆ×ø±êÏµ                                                                     */

	/************************************************************************/
void CKettleDoc::OnOpen3drma() 
{	// TODO: Add your command handler code here
	thd.pointnumber=thd.Read3D_RMA();
 	opensuccess=TRUE;
    UpdateAllViews(NULL);
}

void CKettleDoc::OnSave3drma() 
{
	// TODO: Add your command handler code here
	thd.Save3D_RMAData(thd.pointnumber);
    UpdateAllViews(NULL);
}

void CKettleDoc::OnProcess() 
{
	// TODO: Add your command handler code here
	/*you may process the 3D data here,and redisplay it using OnOpen3drma()*/
}
