#if !defined(AFX_THREED_H__10B94B1D_F58E_4DCD_89C7_6F9637E487E9__INCLUDED_)
#define AFX_THREED_H__10B94B1D_F58E_4DCD_89C7_6F9637E487E9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#define MAXPOINTNUMBER 20000
// ThreeD.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// ThreeD document


class ThreeD : public CDocument
{
//protected:
public:
	ThreeD();           // protected constructor used by dynamic creation
	ThreeD(ThreeD &thd);
	DECLARE_DYNCREATE(ThreeD)

// Attributes
public:
	bool opensuccess;
	short int data[MAXPOINTNUMBER][3];
	float data1[MAXPOINTNUMBER][3];
	int pointnumber;
	long double max;
	float center[3];
    CString m_FileName;
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(ThreeD)
	public:
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:

	BOOL Save3D_RMAData(int pointnumber);
	BOOL Read3D_RMA();
	virtual ~ThreeD();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(ThreeD)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_THREED_H__10B94B1D_F58E_4DCD_89C7_6F9637E487E9__INCLUDED_)
