// KettleView.cpp : implementation of the CKettleView class
//

#include "stdafx.h"
#include "Kettle.h"
#include "math.h"
#include "KettleDoc.h"
#include "KettleView.h"
#include "ThreeD.h"
////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////

// Standard include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
// Windows include files 
#ifdef _WIN32
#include <windows.h>
#endif
// OpenGL include files 
#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glut.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
extern ThreeD thd;
/////////////////////////////////////////////////////////////////////////////
// CKettleView

IMPLEMENT_DYNCREATE(CKettleView, CView)

BEGIN_MESSAGE_MAP(CKettleView, CView)
	//{{AFX_MSG_MAP(CKettleView)
	ON_COMMAND(IDM_ROTATE, OnRotate)
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_WM_TIMER()
	ON_WM_SIZE()
	ON_WM_MOUSEWHEEL()
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CView::OnFilePrintPreview)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CKettleView construction/destruction

CKettleView::CKettleView()
{
	// TODO: add construction code here
	ThreeD thd;
	m_xRotation=0.0f,m_xTranslation=0.0f ;
	m_yRotation=0.0f ,m_yTranslation=0.0f;
	m_zRotation=0.0f ,m_zTranslation=0.0f;
	scale=1.0f;
	m_LeftButtonDown=FALSE, m_RightButtonDown=FALSE;
	m_bRotation=FALSE;
	hDC=NULL; // Private GDI Device Context
	hWnd=NULL; // Holds Our Window Handle
}

CKettleView::~CKettleView()
{
}

BOOL CKettleView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs
    cs.style|=WS_CLIPSIBLINGS|WS_CLIPCHILDREN;
	return CView::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// CKettleView drawing

void CKettleView::OnDraw(CDC* pDC)
{
	CKettleDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	// TODO: add draw code for native data here
	hWnd=GetSafeHwnd();
    hDC=::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLRC);
	if(pDoc->opensuccess){
		DrawScene(pDoc);
	}
	glEnd();
    glPopMatrix();
	glFlush();
	SwapBuffers(hDC);

}

/////////////////////////////////////////////////////////////////////////////
// CKettleView printing

BOOL CKettleView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CKettleView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CKettleView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}
void CKettleView::DrawScene(CKettleDoc* pDoc)
{
    glClearColor(0.9f,0.9f,0.9f,0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear The Screen And The Depth Buffer
	glLoadIdentity(); // Reset The View
	glPushMatrix(); 
    GLfloat light_position[]={0.0f,-1.0f,-1.0f,1.0f};
    glLightfv(GL_LIGHT0,GL_POSITION,light_position);
    GLfloat light_ambient[]={0.5f,0.5f,1.0f,1.0f};
    glLightfv(GL_LIGHT0,GL_AMBIENT,light_ambient);
    GLfloat mat_emission[]={1.0f,1.0f,1.0f,1.0f};
	glMaterialfv(GL_FRONT,GL_EMISSION,mat_emission);
    GLfloat lmodel_ambient[]={0.2f,0.2f,0.2f,1.0f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,lmodel_ambient);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glTranslated( m_xTranslation, m_yTranslation,m_zTranslation);
	glRotatef(m_xRotation,1.0f,0.0f,0.0f);
	glRotatef(m_yRotation,0.0f,1.0f,0.0f);
	glRotatef(m_zRotation,0.0f,0.0f,1.0f);
	glScalef(scale,scale,scale);	
//	glScalef(15.0/pDoc->thd.max,15.0/pDoc->thd.max,15.0/pDoc->thd.max);
//	glTranslated( -pDoc->thd.center[0], -pDoc->thd.center[1],-pDoc->thd.center[2]);


	//short int data[MAXPOINTNUMBER][3]=thd.data;
	//float *data1=&(thd.data1[0][0]);
	glBegin(GL_POINTS);
	glColor3f(0.9f,0.1f,0.2f);
	for(long int i=0;i<pDoc->thd.pointnumber;i++){
		glVertex3f(pDoc->thd.data1[i][0],pDoc->thd.data1[i][1],pDoc->thd.data1[i][2]);
	}
	glEnd();
    glPopMatrix();
	glFlush();
}//在其中添加图形绘制函数。
int CKettleView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	PIXELFORMATDESCRIPTOR pfd=
	{  
		sizeof(PIXELFORMATDESCRIPTOR ),
			1,
			PFD_DRAW_TO_WINDOW|
			 PFD_SUPPORT_OPENGL|
			 PFD_DOUBLEBUFFER,
            PFD_TYPE_RGBA,
			24,
			0,0,0,0,0,0,
			0,
			0,
			0,
			0,0,0,0,
			32,
			0,
			0,
			PFD_MAIN_PLANE,
			0,
			0,0,0
    };
	CClientDC dc(this);
	int pixelformat=ChoosePixelFormat(dc.GetSafeHdc(),&pfd);
	if(SetPixelFormat(dc.GetSafeHdc(),pixelformat,&pfd)==FALSE)
	{
		MessageBox("SetPixelFormat failed");
	    return -1;
	}
    m_hGLRC=wglCreateContext(dc.GetSafeHdc()); 
	
	return 0;
}
//用来说明OpenGL的内部规定。
void CKettleView::OnDestroy()
{
	if(wglGetCurrentContext()!=NULL)
        wglMakeCurrent(NULL,NULL);
    if(m_hGLRC!=NULL)
	{
		wglDeleteContext(m_hGLRC);
	    m_hGLRC =NULL;
	}
	CView::OnDestroy();
}

void CKettleView::OnLButtonDown(UINT nFlags, CPoint point)
{
	m_LeftButtonDown = TRUE;
	mousex=(float)point.x;
	mousey=(float)point.y;	
	CView::OnLButtonDown(nFlags, point);
}
//规定当鼠标左键被按下时进行的操作。
void CKettleView::OnLButtonUp(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_LeftButtonDown=FALSE;	
	CView::OnLButtonUp(nFlags, point);
}
void CKettleView::OnMouseMove(UINT nFlags, CPoint point)
{
	 if(m_LeftButtonDown){ 
        m_xRotation+=(float)(point.y-mousey )/1.0f;
		m_yRotation-=(float)(point.x-mousex )/1.0f;
		mousex=(float)point.x;
	    mousey=(float)point.y;
		InvalidateRect(NULL,FALSE);
	 }
	 if(m_RightButtonDown){
		 m_xTranslation+=(float)(point.x-mousex)/400.0f;
		 m_yTranslation-=(float)(point.y-mousey)/300.0f;
		 mousex=(float)point.x;
		 mousey=(float)point.y;
		 InvalidateRect(NULL,FALSE);
	 }
	CView::OnMouseMove(nFlags, point);
}
//规定当鼠标移动时进行的操作。
void CKettleView::OnRButtonDown(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_RightButtonDown = TRUE;
	mousex=(float)point.x;
	mousey=(float)point.y;  
	
	CView::OnRButtonDown(nFlags, point);
}
void CKettleView::OnRButtonUp(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_RightButtonDown=FALSE;	 
	CView::OnRButtonUp(nFlags, point);
}



void CKettleView::OnTimer(UINT nIDEvent) 
{
	// TODO: Add your message handler code here and/or call default
	switch(nIDEvent)
	{
	case 0:
		break;
    case 1:
		m_xRotation-=0.5f; 
		m_yRotation-=0.5f;
		m_zRotation-=0.5f;
		InvalidateRect(NULL,FALSE);
		break;
    }
	CView::OnTimer(nIDEvent);
}

void CKettleView::OnRotate() 
{
	m_bRotation=!m_bRotation;
	if(m_bRotation)
	 SetTimer(1,10,NULL);
	else
	  KillTimer(1);
}

void CKettleView::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);
		
	if(cy>0)
	{
		CClientDC dc(this);
		wglMakeCurrent(dc.GetSafeHdc(),m_hGLRC);
		glViewport(0,0,cx,cy);
		wglMakeCurrent(NULL,NULL);
	}
}
/////////////////////////////////////////////////////////////////////////////
// CKettleView diagnostics

#ifdef _DEBUG
void CKettleView::AssertValid() const
{
	CView::AssertValid();
}

void CKettleView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CKettleDoc* CKettleView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CKettleDoc)));
	return (CKettleDoc*)m_pDocument;
}
#endif //_DEBUG


BOOL CKettleView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt) 
{
	// TODO: Add your message handler code here and/or call default
	scale*=exp(zDelta/120*0.1);
    InvalidateRect(NULL,FALSE);
	return CView::OnMouseWheel(nFlags, zDelta, pt);
}
