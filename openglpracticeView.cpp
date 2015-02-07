// openglpracticeView.cpp : implementation of the COpenglpracticeView class
//

#include "stdafx.h"
#include "openglpractice.h"
#include "windows.h" // Header File For Windows
#include "openglpracticeDoc.h"
#include "openglpracticeView.h"
#define  StripeImageWidth  32
float   m_xRotation=0.0,m_xTranslation ;
float   m_yRotation=0.0 ,m_yTranslation;
char   m_LeftButtonDown, m_RightButtonDown;
float mousex,mousey;  
float x,y;  
float StripeImage[3*StripeImageWidth];
int j,polymode;
float GenParams[]={1.0,1.0,1.0,0.0}; 


bool keys[256]; // Array Used For The Keyboard Routine
bool active=TRUE; // Window Active Flag
bool fullscreen=TRUE; // Fullscreen Flag Set To TRUE By Default
GLfloat rtri=30.0f; // Angle For The Triangle ( NEW )
GLfloat rquad=30.0f; // Angle For The Quad ( NEW )
GLfloat xrot; // X Rotation ( NEW )
GLfloat yrot; // Y Rotation ( NEW )
GLfloat zrot; // Z Rotation ( NEW )
GLuint texture[3]; // Storage For One Texture ( NEW )
LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM); // Declaration For WndProc


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeView

IMPLEMENT_DYNCREATE(COpenglpracticeView, CView)

BEGIN_MESSAGE_MAP(COpenglpracticeView, CView)
	//{{AFX_MSG_MAP(COpenglpracticeView)
	ON_WM_DESTROY()
	ON_WM_TIMER()
		ON_WM_CREATE()
	ON_WM_SIZE()
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CView::OnFilePrintPreview)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeView construction/destruction

COpenglpracticeView::COpenglpracticeView()
{
	// TODO: add construction code here

}

COpenglpracticeView::~COpenglpracticeView()
{
}

BOOL COpenglpracticeView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeView drawing

void COpenglpracticeView::OnDraw(CDC* pDC)
{
	COpenglpracticeDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	// TODO: add draw code for native data here
	
	hWnd=GetSafeHwnd();
    hDC=::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLRC);
	DrawGLScene();
	wglMakeCurrent(NULL,NULL);
	SwapBuffers(hDC);

}

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeView printing

BOOL COpenglpracticeView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void COpenglpracticeView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void COpenglpracticeView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeView diagnostics

#ifdef _DEBUG
void COpenglpracticeView::AssertValid() const
{
	CView::AssertValid();
}

void COpenglpracticeView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

COpenglpracticeDoc* COpenglpracticeView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(COpenglpracticeDoc)));
	return (COpenglpracticeDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// COpenglpracticeView message handlers

int COpenglpracticeView::DrawGLScene(GLvoid)
{
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear The Screen And The Depth Buffer
glLoadIdentity();// Reset The View

glTranslatef(-0.4f,0.0f,-0.0f); // Move Left 1.5 Units And Into The Screen 6.0
glRotatef(rtri,1.0f,1.0f,1.0f); // Rotate The Triangle On The Y axis ( NEW )

glBegin(GL_TRIANGLES); // Drawing Using Triangles
glColor3f(1.0f,0.0f,0.0f); // Set The Color To Red
glVertex3f( 0.0f, 0.2f, 0.0f); // Top
glColor3f(0.0f,1.0f,0.0f); // Set The Color To Green
glVertex3f(-0.2f,-0.2f, 0.2f); // Bottom Left
glColor3f(0.0f,0.0f,1.0f); // Set The Color To Blue
glVertex3f( 0.2f,-0.2f, 0.2f); // Bottom Right

glColor3f(1.0f,0.0f,0.0f); // Set The Color To Red
glVertex3f( 0.0f, 0.2f, 0.0f); // Top
glColor3f(0.0f,1.0f,0.0f); // Set The Color To Green
glVertex3f(-0.2f,-0.2f, 0.2f); // Bottom Left
glColor3f(0.0f,0.0f,1.0f); // Set The Color To Blue
glVertex3f( -0.0f,-0.2f, -0.2f); // Bottom Right

glColor3f(1.0f,0.0f,0.0f); // Set The Color To Red
glVertex3f(-0.2f,-0.2f, 0.2f); // Top
glColor3f(0.0f,1.0f,0.0f); // Set The Color To Green
glVertex3f(0.2f,-0.2f, 0.2f); // Bottom Left
glColor3f(0.0f,0.0f,1.0f); // Set The Color To Blue
glVertex3f( -0.0f,-0.2f, -0.2f); // Bottom Right

glColor3f(1.0f,0.0f,0.0f); // Set The Color To Red
glVertex3f( 0.0f, 0.2f, 0.0f); // Top
glColor3f(0.0f,1.0f,0.0f); // Set The Color To Green
glVertex3f(0.2f,-0.2f, 0.2f); // Bottom Left
glColor3f(0.0f,0.0f,1.0f); // Set The Color To Blue
glVertex3f( -0.0f,-0.2f, -0.2f); // Bottom Right
glEnd(); // Finished Drawing The Triangle




glLoadIdentity(); // Reset The Current Modelview Matrix because of rotation
glTranslatef(0.5f,0.0f,0.0f); // Move Right 3 Units
//glRotatef(rquad,1.0f,1.0f,1.0f); // Rotate The Quad On The X axis ( NEW )
glRotatef(xrot,1.0f,0.0f,0.0f); // Rotate On The X Axis
glRotatef(yrot,0.0f,1.0f,0.0f); // Rotate On The Y Axis
glRotatef(zrot,0.0f,0.0f,1.0f); // Rotate On The Z Axis
glBindTexture(GL_TEXTURE_2D, texture[2]); // Select Our Texture

glBegin(GL_QUADS); // Draw A Quad
glColor3f(0.5f,0.8f,0.0f);
glTexCoord2f(-1.0f, 1.0f);glVertex3f(-0.2f, 0.2f, 0.2f); // Top Left
glTexCoord2f(1.0f, 1.0f);glVertex3f( 0.2f, 0.2f, 0.2f); // Top Right
glTexCoord2f(1.0f, -1.0f);glVertex3f( 0.2f,-0.2f, 0.2f); // Bottom Right
glTexCoord2f(-1.0f, -1.0f);glVertex3f(-0.2f,-0.2f, 0.2f); // Bottom Left
//glColor3f(0.0f,0.8f,0.2f);

glVertex3f( 0.2f,-0.2f,-0.2f); // Bottom Left Of The Quad (Back)
glVertex3f(-0.2f,-0.2f,-0.2f); // Bottom Right Of The Quad (Back)
glVertex3f(-0.2f, 0.2f,-0.2f); // Top Right Of The Quad (Back)
glVertex3f( 0.2f, 0.2f,-0.2f); // Top Left Of The Quad (Back)
//glColor3f(0.1f,0.1f,0.9f);

glVertex3f( 0.2f,0.2f,0.2f); // Bottom Left Of The Quad (Back)
glVertex3f(0.2f,0.2f,-0.2f); // Bottom Right Of The Quad (Back)
glVertex3f(0.2f, -0.2f,-0.2f); // Top Right Of The Quad (Back)
glVertex3f( 0.2f,-0.2f,0.2f); // Top Left Of The Quad (Back)
//glColor3f(0.9f,0.8f,0.0f);

glVertex3f(-0.2f,-0.2f,0.2f); // Bottom Left Of The Quad (Back)
glVertex3f(-0.2f,-0.2f,-0.2f); // Bottom Right Of The Quad (Back)
glVertex3f(-0.2f, 0.2f,-0.2f); // Top Right Of The Quad (Back)
glVertex3f(-0.2f, 0.2f,0.2f); // Top Left Of The Quad (Back)
//glColor3f(0.9f,0.0f,0.1f);

glVertex3f( 0.2f,0.2f,-0.2f); // Bottom Left Of The Quad (Back)
glVertex3f(-0.2f,0.2f,0.2f); // Bottom Right Of The Quad (Back)
glVertex3f(-0.2f,0.2f,-0.2f); // Top Right Of The Quad (Back)
glVertex3f( 0.2f,0.2f,0.2f); // Top Left Of The Quad (Back)
//glColor3f(0.5f,0.0f,0.9f);

glVertex3f( 0.2f,-0.2f,-0.2f); // Bottom Left Of The Quad (Back)
glVertex3f(-0.2f,-0.2f,-0.2f); // Bottom Right Of The Quad (Back)
glVertex3f(-0.2f,-0.2f,0.2f); // Top Right Of The Quad (Back)
glVertex3f( 0.2f,-0.2f,0.2f); // Top Left Of The Quad (Back)
glEnd(); // Done Drawing The Quad



 SetTimer(1,10,NULL);	
return 1;//always running
}

int COpenglpracticeView::OnCreate(LPCREATESTRUCT lpCreateStruct)
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

void COpenglpracticeView::OnDestroy()
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
void COpenglpracticeView::OnTimer(UINT nIDEvent) 
{
	// TODO: Add your message handler code here and/or call default
	switch(nIDEvent)
	{
	case 0:
		break;
    case 1:
		xrot+=3.0f; // X Axis Rotation
		yrot+=2.0f; // Y Axis Rotation
		zrot+=4.0f; // Z Axis Rotation
		rtri+=2.0f; // Increase The Rotation Variable For The Triangle ( NEW )
		rquad-=1.5f; // Decrease The Rotation Variable For The Quad ( NEW )
		InvalidateRect(NULL,FALSE);
		break;
     
    }

	CView::OnTimer(nIDEvent);
}

GLvoid COpenglpracticeView::OnSize(UINT nType, GLsizei width, GLsizei height) 
{
	if(height>0)
	{
		CClientDC dc(this);
		wglMakeCurrent(dc.GetSafeHdc(),m_hGLRC);
		glViewport(0,0,width, height);
		wglMakeCurrent(NULL,NULL);
	}
	// TODO: Add your message handler code here
}
