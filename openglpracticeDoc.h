// openglpracticeDoc.h : interface of the COpenglpracticeDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_OPENGLPRACTICEDOC_H__0B4324E2_2AFE_4FFB_9111_C2B76D164F2A__INCLUDED_)
#define AFX_OPENGLPRACTICEDOC_H__0B4324E2_2AFE_4FFB_9111_C2B76D164F2A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


class COpenglpracticeDoc : public CDocument
{
protected: // create from serialization only
	COpenglpracticeDoc();
	DECLARE_DYNCREATE(COpenglpracticeDoc)

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(COpenglpracticeDoc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~COpenglpracticeDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(COpenglpracticeDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OPENGLPRACTICEDOC_H__0B4324E2_2AFE_4FFB_9111_C2B76D164F2A__INCLUDED_)
