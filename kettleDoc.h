 // kettleDoc.h : interface of the CKettleDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_KETTLEDOC_H__96235AA8_A39B_4E3A_8DE7_3FCA1236D4F5__INCLUDED_)
#define AFX_KETTLEDOC_H__96235AA8_A39B_4E3A_8DE7_3FCA1236D4F5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "ThreeD.h"


class CKettleDoc : public CDocument
{
protected: // create from serialization only
	CKettleDoc();
	DECLARE_DYNCREATE(CKettleDoc)

// Attributes
public:
bool opensuccess;
ThreeD thd;
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CKettleDoc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CKettleDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

public:
	afx_msg void OnOpen3drma();
// Generated message map functions
protected:
	//{{AFX_MSG(CKettleDoc)
	afx_msg void OnSave3drma();
	afx_msg void OnProcess();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_KETTLEDOC_H__96235AA8_A39B_4E3A_8DE7_3FCA1236D4F5__INCLUDED_)
