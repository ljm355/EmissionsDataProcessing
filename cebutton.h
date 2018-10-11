#pragma once
//#include "../celengine/texture.h"
#include "cetexture.h"
#include "gl/glew.h"
#include <vector>
#include <osg/Vec4>
#include <osg/Vec2>
#include <osg/Vec3>
#include <osg/Matrix>
#include <map>

class CECallBackBase : public osg::Referenced
{
public:
	CECallBackBase(void){}
	virtual ~CECallBackBase(void){}
};

class CECallBackManger : public osg::Referenced
{

public:
	void addCallBack(CECallBackBase* callback);
	void removeCallBack(CECallBackBase* callback);
	typedef std::map<CECallBackBase*,CECallBackBase*> Map;
protected:
	std::map<CECallBackBase*,CECallBackBase*> g_mCallBackMap;
};

struct VertexUVXYZ
{
	float tu, tv;
	float x, y, z;
	VertexUVXYZ()
	{

	}
	VertexUVXYZ(float u,float v,float px,float py,float pz)
		: tu(u),tv(v),x(px),y(py),z(pz)
	{

	}
};


struct CERect
{
	int x,y,w,h;
	CERect()
	{
		x = 0;y=0;w=0;h=0;
	}
	CERect(int px,int py,int width,int height):
	x(px),y(py),w(width),h(height)
	{

	}
}; 
struct CERectF
{
	float x,y,w,h;
	CERectF()
	{
		x = 0;y=0;w=0;h=0;
	}
	CERectF(float px,float py,float width,float height):
	x(px),y(py),w(width),h(height)
	{

	}
}; 
struct CEButtonState
{
	osg::Vec4 highlightColor;
	osg::Vec4 textureRect;
    osg::Matrix transform;
	std::string name;
	int id;
	CEButtonState()
	{
		textureRect = osg::Vec4(0,0,1,1);
		transform = osg::Matrix::identity();
		highlightColor = osg::Vec4(0,0,0,0);
		name = "";
	}
};

class CEButtonCallBack : public CECallBackBase, public osg::Referenced
{
public:
	virtual void onMouseHover(int x,int y,int button,void* object){}
	virtual void onMouseEnter(int x,int y,int button,void* object){}
	virtual void onMouseLeave(int x,int y,int button,void* object){}
	virtual void onMouseDown(int x,int y,int button,void* object){}
	virtual void onMouseUp(int x,int y,int button,void* object){}
};



class CEButton : public osg::Referenced
{
public:
	class CEButtonVisitor : public osg::Referenced
	{
     public:
		virtual void applyButton(CEButton* button) = 0;
	};
	CEButton(void);
	~CEButton(void);
	virtual void traverse(CEButtonVisitor* visitor);
	virtual void getNormalizedRect(VertexUVXYZ* pVertex,const int& screenWidth,const int&  screenHeight);
	virtual void draw(const int& screenWidth,const int& screenHeight);
	virtual void drawThis(const int& screenWidth,const int&  screenHeight);
	virtual void drawChildren(const int& screenWidth,const int&  screenHeight);
	virtual CERect Rect() const { return g_mRect; }
	//virtual void addCallBack(CECallBackBase* callback );
	virtual osg::Vec2 screenPos();
	virtual void resize(int x,int y,int w,int h);
	virtual void addChild(CEButton* child){g_mChildren.push_back(child);child->parent(this);}
	virtual int numChildren(){return g_mChildren.size();}
	virtual CEButton* getChild(int idx){return g_mChildren[idx];}
	virtual void getAbsolutePos(int& x,int& y);
	virtual CEButton* parent() const { return g_pParent; }
	virtual void parent(CEButton* val) { g_pParent = val; }
	virtual CEButtonState& defaultState() { return g_mDefaultState; }
	virtual CEButtonState& currentState() { return g_mCurrentState; }
	virtual void setDefaultState(const CEButtonState& val) { g_mDefaultState = val; }
	virtual void applyDefaultState(){g_mCurrentState = g_mDefaultState;}
	virtual void setCurrentState(const CEButtonState& state){g_mCurrentState=state;}
	virtual void addState(const CEButtonState& state){g_mStates.push_back(state);}
	virtual CEButtonState& getState(int idx){return g_mStates[idx];}
	virtual int numStates(const CEButtonState& state){return g_mStates.size();}
	virtual void applyState(int idx){g_mCurrentState=g_mStates[idx];}
	virtual CEButton* onMouseDown(int mousex, int mousey,int button);
	virtual CEButton* onMouseUp(int mousex, int mousey,int button);
	virtual CEButton* onMouseMove(int mousex, int mousey,int button,CECallBackManger::Map& callbacks);
	virtual bool hasMouseEntered() const { return g_mHasMouseEntered; }
	virtual CETexture* getTexture() const { return g_pTexture; }
	virtual void setTexture(CETexture* val);
	bool isOn() const { return g_mIsOn; }
	void isOn(bool val) { g_mIsOn = val; }
	float alpha() const { return g_mAlpha; }
	void alpha(float val) { g_mAlpha = val; }
	void setShaderFromString(std::string fragmentShader,std::string vertexShader=defaultVertexShader());
	static std::string readTextFromFile(const char *fn);
	static GLuint defaultShaderProgam();
	std::string name() const { return g_mName; }
	void name(std::string val) { g_mName = val; }
	bool IsAutoFill() const { return g_mIsAutoFill; }
	void IsAutoFill(bool val) { g_mIsAutoFill = val; }
protected:
	static std::string defaultVertexShader();
	static std::string defaultFragmentShader();
	void freeShader();
	void setShaderFromFile(std::string vertexShaderFile,std::string fragmentShaderFile="");
	virtual bool isMouseOn(int mousex, int mousey);
protected: 
	osg::ref_ptr<CETexture> g_pTexture;
	CERect g_mRect;
	float g_mAlpha;
	CEButton* g_pParent;
	std::vector<CEButton*> g_mChildren;
	CEButtonState g_mDefaultState;
	CEButtonState g_mCurrentState;
	std::vector<CEButtonState> g_mStates;
	bool g_mHasMouseEntered;
	bool g_mIsOn;
	GLuint g_pProgramHandle;
	GLuint g_pVShader;
	GLuint g_pFShader;
	GLuint g_pDefaultProgramHandle;
	std::string g_mName;
	bool g_mIsAutoFill;
	CERectF g_mNormalizedRect;

};

class CEButtonManager : public CEButton,public CECallBackManger,  public osg::Referenced
{
public:
	virtual CEButton* onMouseDown(int mousex, int mousey,int button);
	virtual CEButton* onMouseUp(int mousex, int mousey,int button);
	virtual CEButton* onMouseMove(int mousex, int mousey,int button);
	virtual void addCallBack(CECallBackBase* callback );
};

class CEButtonFading : public CEButton::CEButtonVisitor
{
public:
	enum FadeState
	{
		FadingIn,
		FadingOut,
		Finished
	};

	virtual float alpha() const { return g_mAlpha; }
	virtual FadeState state() const { return g_mState; }
	virtual float alphaStep() const { return g_mAlphaStep; }
	virtual void alphaStep(float val) { g_mAlphaStep = val; }
	virtual void fadeIn();
	virtual void fadeOut();
	virtual void fade();
	//void Alpha(float val) { g_mAlpha = val; }
	virtual void applyButton(CEButton* button);
	virtual void attachMenuPage(CEButtonManager* menuPage);
	CEButtonFading(CEButtonManager* menuPage);
private:
	FadeState g_mState;
	float g_mAlpha;
	CEButtonManager* g_pMenuPage;
	float g_mAlphaStep;

};
