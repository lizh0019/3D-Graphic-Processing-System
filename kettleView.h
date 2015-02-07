// kettleView.h : interface of the CKettleView class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_KETTLEVIEW_H__B652AC0D_F35B_4FBE_8433_6A153BF20D39__INCLUDED_)
#define AFX_KETTLEVIEW_H__B652AC0D_F35B_4FBE_8433_6A153BF20D39__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000



class CKettleView : public CView
{
protected: // create from serialization only
	CKettleView();
	DECLARE_DYNCREATE(CKettleView)

// Attributes
public:
	CKettleDoc* GetDocument();
	volatile float   m_xRotation,m_xTranslation;
	volatile float   m_yRotation,m_yTranslation;
	volatile float   m_zRotation,m_zTranslation;
	volatile bool   m_LeftButtonDown, m_RightButtonDown;
	volatile float mousex,mousey;  
	volatile float scale;
	volatile bool m_bRotation;
	HDC hDC; // Private GDI Device Context
	HWND hWnd; // Holds Our Window Handle
// Operations
public:
	// GLUT variables 
	HGLRC m_hGLRC;
	void DrawScene(CKettleDoc* pDoc);
	void OnRotate();
// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CKettleView)
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
	virtual ~CKettleView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	int OnCreate(LPCREATESTRUCT lpCreateStruct);//用来说明OpenGL的内部规定。
	void OnDestroy();
// Generated message map functions
protected:
	//{{AFX_MSG(CKettleView)
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);//规定当鼠标左右键被按下时进行的操作。
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);//规定当鼠标移动时进行的操作。
	afx_msg void OnTimer(UINT nIDEvent);
	afx_msg void OnSize(UINT nType,int cx,int cy);	//规定视区的大小和投影变换。
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in kettleView.cpp
inline CKettleDoc* CKettleView::GetDocument()
   { return (CKettleDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_KETTLEVIEW_H__B652AC0D_F35B_4FBE_8433_6A153BF20D39__INCLUDED_)
