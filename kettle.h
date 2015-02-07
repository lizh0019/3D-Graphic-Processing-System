// kettle.h : main header file for the KETTLE application
//

#if !defined(AFX_KETTLE_H__0E29EED2_B40E_4D73_978F_D4EA15A9D88D__INCLUDED_)
#define AFX_KETTLE_H__0E29EED2_B40E_4D73_978F_D4EA15A9D88D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// CKettleApp:
// See kettle.cpp for the implementation of this class
//

class CKettleApp : public CWinApp
{
public:
	CKettleApp();


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CKettleApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation
	//{{AFX_MSG(CKettleApp)
	afx_msg void OnAppAbout();
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_KETTLE_H__0E29EED2_B40E_4D73_978F_D4EA15A9D88D__INCLUDED_)
