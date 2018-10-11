#include "cebutton.h"
#include "glut.h"

CEButton::CEButton(void)
	:g_pParent(NULL),g_mHasMouseEntered(false),g_pTexture(NULL),g_mAlpha(1)
{
	g_pFShader = 0;
	g_pVShader = 0;
	g_pProgramHandle = 0;
	g_mDefaultState.transform = osg::Matrix::identity();
	applyDefaultState();
	g_mIsOn = true;
	g_mIsAutoFill = false;
}

CEButton::~CEButton(void)
{
	for (int i=0;i<g_mChildren.size();i++)
	{
		if(g_mChildren[i])
			delete g_mChildren[i];
	}
	g_mChildren.clear();
	freeShader();
}

void CEButton::setTexture( CETexture* val )
{
	g_pTexture = val;

}

void CEButton::drawThis(const int& screenWidth,const int&  screenHeight)
{

	VertexUVXYZ pVertex[4];
	getNormalizedRect(pVertex,screenWidth,screenHeight);
	//glDisable(GL_CULL_FACE);
	glEnable(GL_TEXTURE_2D);				// Depth Buffer Setup
	glDisable(GL_DEPTH_TEST);	
	glEnable(GL_BLEND) ;
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();										// Reset The Modelview Matrix

	if(g_pTexture)
		g_pTexture->bind();
	    
	glInterleavedArrays( GL_T2F_V3F, 0, pVertex );
	glActiveTexture(GL_TEXTURE0);

	GLuint shaderProgram = g_pProgramHandle;
	if(!shaderProgram)
		shaderProgram = defaultShaderProgam();

	glUseProgram(shaderProgram);
	int loc = glGetUniformLocation(shaderProgram, "texture0");
	if(loc >-1)
		glUniform1i(loc, 0);
	loc = glGetUniformLocation(shaderProgram, "texrect");
	if(loc >-1)
		glUniform4f(loc, g_mCurrentState.textureRect.x(), g_mCurrentState.textureRect.y(),
		g_mCurrentState.textureRect.z(), g_mCurrentState.textureRect.w());
	loc = glGetUniformLocation(shaderProgram, "color");
	if(loc >-1)
		glUniform4f(loc, g_mCurrentState.highlightColor.x(), g_mCurrentState.highlightColor.y(),
		g_mCurrentState.highlightColor.z(), g_mCurrentState.highlightColor.w());
	loc = glGetUniformLocation(shaderProgram, "alpha");
	if(loc >-1)
		glUniform1f(loc, g_mAlpha);

	glDrawArrays( GL_QUADS, 0, 4 );
	glUseProgram(0);


	glEnable(GL_DEPTH_TEST);	
	glDisable(GL_BLEND) ;
	glDisable(GL_VERTEX_ARRAY);

}

void CEButton::drawChildren( const int& screenWidth,const int& screenHeight )
{
	for (int i=0;i<g_mChildren.size();i++)
	{
    
		//if(g_mChildren[i])
		g_mChildren[i]->draw(screenWidth,screenHeight);
	}
}

void CEButton::draw(const int& screenWidth,const int& screenHeight)
{
	if(!g_mIsOn)
		return;
	drawThis(screenWidth,screenHeight);	
	drawChildren(screenWidth,screenHeight);
}

CEButton* CEButton::onMouseDown( int mousex, int mousey,int button )
{  
	for (int i=0;i<g_mChildren.size();i++)
	{
		CEButton* sel = g_mChildren[i]->onMouseDown(mousex,mousey,button);
		if(sel)
			return sel;
		//if(g_mChildren[i])
		//	g_mChildren[i]->draw(screenWidth,screenHeight);
	}
	if(isMouseOn(mousex,mousey))
	{  
		return this;  
	}  
	return NULL;  
}  

CEButton* CEButton::onMouseMove( int mousex, int mousey,int button,CECallBackManger::Map& callbacks )
{

	for (int i=0;i<g_mChildren.size();i++)
	{
		CEButton* sel = g_mChildren[i]->onMouseMove(mousex,mousey,button,callbacks);
		//if(sel)
		//	return sel;
	}
	//printf("move\n");
	if(isMouseOn(mousex,mousey))
	{  
		if(!g_mHasMouseEntered)
		{
			CECallBackManger::Map::iterator iter = callbacks.begin();
			while(iter != callbacks.end())
			{
				if(dynamic_cast<CEButtonCallBack*>(iter->second))
					(dynamic_cast<CEButtonCallBack*>(iter->second))->onMouseEnter(mousex,mousey,button,this);
				iter++;
			}
		}
		else
		{
			CECallBackManger::Map::iterator iter = callbacks.begin();
			while(iter != callbacks.end())
			{
				if(dynamic_cast<CEButtonCallBack*>(iter->second))
					(dynamic_cast<CEButtonCallBack*>(iter->second))->onMouseHover(mousex,mousey,button,this);
				iter++;
			}
		}
		g_mHasMouseEntered = true;
		return this;  
	}  
	else
	{
		if(g_mHasMouseEntered)
		{
			CECallBackManger::Map::iterator iter = callbacks.begin();
			while(iter != callbacks.end())
			{
				if(dynamic_cast<CEButtonCallBack*>(iter->second))
					(dynamic_cast<CEButtonCallBack*>(iter->second))->onMouseLeave(mousex,mousey,button,this);
				iter++;
			}
		}
		g_mHasMouseEntered = false;
	}
	return NULL;  
}

CEButton* CEButton::onMouseUp(int mousex, int mousey,int button)
{  

	for (int i=0;i<g_mChildren.size();i++)
	{
		CEButton* sel = g_mChildren[i]->onMouseUp(mousex,mousey,button);
		if(sel)
			return sel;
		//if(g_mChildren[i])
		//	g_mChildren[i]->draw(screenWidth,screenHeight);
	}
	if(isMouseOn(mousex,mousey))
	{  
		return this;  
	}  
	return NULL;  
}  

void CEButton::resize( int x,int y,int w,int h )
{  
	CERect oldRect = g_mRect;
	if(oldRect.w == 0 && oldRect.h == 0)
	{
		oldRect = CERect(x,y,w,h);
	}
	g_mRect = CERect(x,y,w,h);
	float wscale = (float)w/oldRect.w;
	float hscale = (float)h/oldRect.h;

	
		for (int i=0;i<g_mChildren.size();i++)
		{
			CEButton& btn = *g_mChildren[i];
			if(btn.IsAutoFill())
			{
				CERect rect = CERect((int)(btn.Rect().x*wscale),(int)(btn.Rect().y*hscale),(int)(btn.Rect().w*wscale),(int)(btn.Rect().h*hscale));
				g_mChildren[i]->resize(rect.x,rect.y,rect.w,rect.h);
			}

		}
	
}

void CEButton::getNormalizedRect( VertexUVXYZ* pVertex,const int& screenWidth,const int&  screenHeight)
{
	int x,y;
	getAbsolutePos(x,y);
	osg::Matrix transform = g_mCurrentState.transform;
	 
	osg::Vec3 origin = osg::Vec3(x+g_mRect.w*0.5,g_mRect.h*0.5+y,0);
	transform = osg::Matrix::translate(-origin) * transform * osg::Matrix::translate(origin);

	std::vector<osg::Vec4> posArr;
	posArr.push_back(osg::Vec4(x,          y,          1,0));
	posArr.push_back(osg::Vec4(x+g_mRect.w,y,          1,0));
	posArr.push_back(osg::Vec4(x+g_mRect.w,y+g_mRect.h,1,0));
	posArr.push_back(osg::Vec4(x,          y+g_mRect.h,1,0));

	pVertex[0]= VertexUVXYZ(0.0f,1.0f, 1.0f,1.0f, 0.0f);
	pVertex[1]= VertexUVXYZ(1.0f,1.0f, 1.0f,1.0f, 0.0f);
	pVertex[2]= VertexUVXYZ(1.0f,0.0f, 1.0f,1.0f, 0.0f);
	pVertex[3]= VertexUVXYZ(0.0f,0.0f, 1.0f,1.0f, 0.0f);

	for (int i=0;i<4;i++)
	{
		osg::Vec3 pos = osg::Vec3(posArr[i].x(),posArr[i].y(),posArr[i].z());
		pos = pos * transform;
		pos.x() = (float)pos.x()/screenWidth;
		pos.y() = (float)pos.y()/screenHeight;

		pos.x() = pos.x() * 2.0 - 1.0;
		pos.y() = pos.y() * 2.0 - 1.0;

		pVertex[i].x = pos.x();
		pVertex[i].y = pos.y();

	}




}

void CEButton::getAbsolutePos( int& x,int& y )
{
	CEButton* parent = g_pParent;
	x = g_mRect.x;
	y = g_mRect.y;
	while(parent)
	{
		x+=parent->Rect().x;
		y+=parent->Rect().y;
		if(parent->parent())
			parent = parent->parent();
		else
			break;
	}
}

osg::Vec2 CEButton::screenPos()
{
	int x,y;
	getAbsolutePos(x,y);
	return osg::Vec2(x,y); 
}

bool CEButton::isMouseOn(int mousex, int mousey)
{
	int x,y;
	getAbsolutePos(x,y);
	int w = glutGet(GLUT_WINDOW_WIDTH);
	int h = glutGet(GLUT_WINDOW_HEIGHT); 
	mousey = h-mousey;  
	if( mousex > x && mousex < x+g_mRect.w &&  
		mousey > y && mousey < y+g_mRect.h )  
	{  
	
		return true;  
	}  
	return false;
}

void CEButton::freeShader()
{
	if(g_pProgramHandle)
	{
		if(g_pFShader)
		{
			glDetachShader(g_pProgramHandle, g_pFShader);
			glDeleteShader(g_pFShader);
		}
		if(g_pVShader)
		{
			glDetachShader(g_pProgramHandle, g_pVShader);
			glDeleteShader(g_pVShader);
		}

	}

	g_pFShader = 0;
	g_pVShader = 0;
	g_pProgramHandle = 0;
}

std::string CEButton::readTextFromFile(const char *fn) {

	FILE *fp;
	char *content = NULL;

	int count=0;

	if (fn != NULL) {
		fp = fopen(fn,"rt");

		if (fp != NULL) {

			fseek(fp, 0, SEEK_END);
			count = ftell(fp);
			rewind(fp);

			if (count > 0) {
				content = (char *)malloc(sizeof(char) * (count+1));
				count = fread(content,sizeof(char),count,fp);
				content[count] = '\0';
			}
			fclose(fp);
		}
	}
	return content;
}

GLuint CEButton::defaultShaderProgam()
{
	static GLuint programHandle = 0;
	static GLuint vsHandle = 0;
	static GLuint psHandle = 0;
	if(programHandle)
		return programHandle;
	vsHandle = glCreateShader( GL_VERTEX_SHADER );
	std::string vsCode = defaultVertexShader();
	std::string fsCode = defaultFragmentShader();
	GLchar* vShaderCode = const_cast<char*>(vsCode.data());
	const GLchar* vCodeArray[1] = {vShaderCode};
	glShaderSource(vsHandle, 1, vCodeArray, NULL);
	glCompileShader(vsHandle);
	//free(vShaderCode);
	GLint logLen;
	glGetShaderiv(vsHandle, GL_INFO_LOG_LENGTH, &logLen);
	if(logLen > 0)
	{
		char *log = (char *)malloc(logLen);
		GLsizei written;
		glGetShaderInfoLog(vsHandle, logLen, &written, log);
		printf("Shader compile error log: %s\n",log);
		free(log);
	}

	psHandle = glCreateShader( GL_FRAGMENT_SHADER );
	//printf("Here\n");


	GLchar* fShaderCode = const_cast<char*>(fsCode.data());;

	const GLchar* fCodeArray[1] = {fShaderCode};

	glShaderSource(psHandle, 1, fCodeArray, NULL);


	glCompileShader(psHandle);

	programHandle = glCreateProgram();

	glAttachShader(programHandle, vsHandle);

	glAttachShader(programHandle, psHandle);
	glLinkProgram(programHandle);
	return programHandle;

}

void CEButton::setShaderFromString(std::string fragmentShader,std::string vertexShader)
{
	freeShader();
	if(vertexShader != "")
	{
		g_pVShader = glCreateShader( GL_VERTEX_SHADER );
		GLchar* vShaderCode = const_cast<char*>(vertexShader.data());
		const GLchar* vCodeArray[1] = {vShaderCode};
		glShaderSource(g_pVShader, 1, vCodeArray, NULL);
		glCompileShader(g_pVShader);
		//free(vShaderCode);
		GLint logLen;
		glGetShaderiv(g_pVShader, GL_INFO_LOG_LENGTH, &logLen);
		if(logLen > 0)
		{
			char *log = (char *)malloc(logLen);
			GLsizei written;
			glGetShaderInfoLog(g_pVShader, logLen, &written, log);
			printf("Shader compile error log: %s\n",log);
			free(log);
		}

	}

	if(fragmentShader != "")
	{
		g_pFShader = glCreateShader( GL_FRAGMENT_SHADER );
		//printf("Here\n");


		GLchar* fShaderCode = const_cast<char*>(fragmentShader.data());;

		const GLchar* fCodeArray[1] = {fShaderCode};

		glShaderSource(g_pFShader, 1, fCodeArray, NULL);


		glCompileShader(g_pFShader);

		//free(fShaderCode);
		//const GLchar* codeArray[] = {shaderCode};
		//Check the compile result
	}

	if(g_pFShader || g_pVShader)
	{
		g_pProgramHandle = glCreateProgram();
		if(0 == g_pProgramHandle)
		{
			fprintf(stderr, "Error creating programHandle.\n");
		}
		if(g_pVShader)
			glAttachShader(g_pProgramHandle, g_pVShader);
		if(g_pFShader)
			glAttachShader(g_pProgramHandle, g_pFShader);
		glLinkProgram(g_pProgramHandle);
		//int loc = glGetUniformLocation(g_pFShader, "texture0");
		//loc = glGetUniformLocation(g_pFShader, "texrect");
	}
	//glUseProgram(programHandle);
}

void CEButton::setShaderFromFile(std::string fragmentShaderFile,std::string vertexShaderFile)
{
	std::string vShaderCode = "";
	std::string fShaderCode = "";
	if(vertexShaderFile != "")
	{
		vShaderCode = readTextFromFile(vertexShaderFile.data());
	}
	if(fragmentShaderFile != "")
	{
		fShaderCode = readTextFromFile(fragmentShaderFile.data());
	}
	setShaderFromString(fShaderCode, vShaderCode);
}

std::string CEButton::defaultVertexShader()
{
	static std::string defaultVS = 
		"uniform vec4 texrect;\n"
		"void main(void)\n"
		"{\n"
		"   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
		"   gl_TexCoord[1] = gl_MultiTexCoord1;\n"
		"   gl_TexCoord[2] = gl_MultiTexCoord2;\n"
		"   gl_TexCoord[3] = gl_MultiTexCoord3;\n"
		"   vec4 texpos = gl_MultiTexCoord0;\n"
		"   texpos.x = texrect.x+texrect.z*texpos.x;\n"
		"   texpos.y = texrect.y+texrect.w*texpos.y;\n"
		"   gl_TexCoord[0] = texpos;\n"  
		"   gl_Position   = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
		"}\n";
	return defaultVS;
}

std::string CEButton::defaultFragmentShader()
{ 
static std::string  defaultFS = 
	"uniform sampler2D texture0;\n"
	"uniform vec4 color;\n"
	"uniform float alpha;\n"
	"uniform vec4 texrect;\n"
	"\n"
	"void main(void) \n"
	"{\n"
	"    gl_FragColor = texture2D( texture0, gl_TexCoord[0].xy)+color; \n"
	"    gl_FragColor.a = gl_FragColor.a * alpha; \n"
	"}\n";
return defaultFS;
}

void CEButton::traverse( CEButtonVisitor* visitor )
{
    visitor->applyButton(this);
	for (int i=0;i<g_mChildren.size();i++)
	{
		if(g_mChildren[i])
			g_mChildren[i]->traverse(visitor);
	}
}

void CEButtonManager::addCallBack( CECallBackBase* callback )
{
	if(!dynamic_cast<CEButtonCallBack*>(callback))
		return;
	//dynamic_cast<ViewerCallBack*>(callback)->setViewer(this);
	CECallBackManger::addCallBack(callback);
}

CEButton* CEButtonManager::onMouseDown( int mousex, int mousey,int button )
{
	if(!g_mIsOn)
		return NULL;
	CEButton* sel = CEButton::onMouseDown(mousex,mousey,button);
	if(sel)
	{
		CECallBackManger::Map::iterator iter = g_mCallBackMap.begin();
		while(iter != g_mCallBackMap.end())
		{
			if(dynamic_cast<CEButtonCallBack*>(iter->second))
				(dynamic_cast<CEButtonCallBack*>(iter->second))->onMouseDown(mousex,mousey,button,sel);
			iter++;
		}
	}
	return sel;
}

CEButton* CEButtonManager::onMouseUp( int mousex, int mousey,int button )
{
	if(!g_mIsOn)
		return NULL;
	CEButton* sel = CEButton::onMouseUp(mousex,mousey,button);
	if(sel)
	{
		CECallBackManger::Map::iterator iter = g_mCallBackMap.begin();
		while(iter != g_mCallBackMap.end())
		{
			if(dynamic_cast<CEButtonCallBack*>(iter->second))
				(dynamic_cast<CEButtonCallBack*>(iter->second))->onMouseUp(mousex,mousey,button,sel);
			iter++;
		}
	}
	return sel;
}

CEButton* CEButtonManager::onMouseMove( int mousex, int mousey,int button )
{
	if(!g_mIsOn)
		return NULL;
	CEButton* sel = CEButton::onMouseMove(mousex,mousey,button,g_mCallBackMap);
	if(sel)
	{
		//CECallBackManger::Map::iterator iter = g_mCallBackMap.begin();
		//while(iter != g_mCallBackMap.end())
		//{
		//	if(dynamic_cast<CEButtonCallBack*>(iter->second))
		//		(dynamic_cast<CEButtonCallBack*>(iter->second))->onMouseHover(mousex,mousey,button,sel);
		//	iter++;
		//}
	}
	return sel;
}

void CECallBackManger::addCallBack( CECallBackBase* callback )
{
	if(!callback)
		return;
	g_mCallBackMap[callback]=callback;
}

void CECallBackManger::removeCallBack( CECallBackBase* callback )
{
	if(!callback)
		return;
	if(g_mCallBackMap.find(callback) != g_mCallBackMap.end())
		g_mCallBackMap.erase(callback);
}

void CEButtonFading::fadeIn()
{
	g_mState = FadingIn;
	g_mAlpha = 0;
}

void CEButtonFading::fadeOut()
{
	g_mState = FadingOut;
	g_mAlpha = 1;
}

void CEButtonFading::fade()
{
	if(g_mState == FadingIn)
	{
		g_pMenuPage->traverse(this);
		g_mAlpha+=g_mAlphaStep;
	}
	else if(g_mState == FadingOut)
	{
		g_pMenuPage->traverse(this);
		g_mAlpha-=g_mAlphaStep;
	}
	if(g_mAlpha>1)
	{
		g_mState = Finished;
		g_pMenuPage->isOn(true);
		g_mAlpha = 1;
		g_pMenuPage->traverse(this);
	}
	if(g_mAlpha < 0)
	{
		g_pMenuPage->isOn(false);
		g_mState = Finished;
		g_mAlpha = 0;
		g_pMenuPage->traverse(this);
	}
}

void CEButtonFading::applyButton( CEButton* button )
{
	button->alpha(g_mAlpha);
}

void CEButtonFading::attachMenuPage( CEButtonManager* menuPage )
{
	g_pMenuPage = menuPage;
}

CEButtonFading::CEButtonFading( CEButtonManager* menuPage )
{
	g_pMenuPage = menuPage;
	g_mAlpha = 1;
	g_mAlphaStep = 0.01;
	g_mState = Finished;
}
