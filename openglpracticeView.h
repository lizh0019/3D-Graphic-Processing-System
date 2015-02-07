// openglpracticeView.h : interface of the COpenglpracticeView class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_OPENGLPRACTICEVIEW_H__1A57A207_9767_402B_B188_7B806FCE4B07__INCLUDED_)
#define AFX_OPENGLPRACTICEVIEW_H__1A57A207_9767_402B_B188_7B806FCE4B07__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

HDC hDC=NULL; // Private GDI Device Context
HGLRC m_hGLRC=NULL; // Permanent Rendering Context
HWND hWnd=NULL; // Holds Our Window Handle
HINSTANCE hInstance; // Holds The Instance Of The Application
class COpenglpracticeView : public CView
{
protected: // create from serialization only
	COpenglpracticeView();
	DECLARE_DYNCREATE(COpenglpracticeView)

// Attributes
public:
	COpenglpracticeDoc* GetDocument();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(COpenglpracticeView)
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
	//}}AFX_VIRTUAL

// Implementation
public:
	int OnCreate(LPCREATESTRUCT lpCreateStruct);
	int DrawGLScene(GLvoid);
	virtual ~COpenglpracticeView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(COpenglpracticeView)
	afx_msg void OnDestroy();
	afx_msg void OnTimer(UINT nIDEvent);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in openglpracticeView.cpp
inline COpenglpracticeDoc* COpenglpracticeView::GetDocument()
   { return (COpenglpracticeDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OPENGLPRACTICEVIEW_H__1A57A207_9767_402B_B188_7B806FCE4B07__INCLUDED_)
