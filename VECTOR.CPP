#include "vector.h"


	extern Vector3 vector;

    Vector3::Vector3(void):x(0),y(0),z(0){}
	Vector3::Vector3(float xx,float yy,float zz){x=xx;y=yy;z=zz;}
	Vector3& Vector3::operator =(Vector3& vector1){
		 x=vector1.x;y=vector1.y;z=vector1.z;
		 return *this;
	}
	Vector3::Vector3(Vector4 &p){x=p.x;y=p.y;z=p.z;}
	Vector3 Vector3::operator + (Vector3& vector1){
		Vector3	vector = *this;
		vector.x=x+vector1.x;
		vector.y=y+vector1.y;
		vector.z=z+vector1.z;
		return	vector;
	}
	Vector3 Vector3::operator - (Vector3& vector1){
		Vector3	vector = *this;
		vector.x=x-vector1.x;
		vector.y=y-vector1.y;
		vector.z=z-vector1.z;
		return	vector;
	}
	float Vector3::operator / (Vector3& vector1){
		Vector3	vector = *this;
		vector.x=x*vector1.x;
		vector.y=y*vector1.y;
		vector.z=z*vector1.z;
		float sum=vector.x+vector.y+vector.z;
		return	sum;
	}
	Vector3 Vector3::operator / (float dValue){
		Vector3 vector = *this;
		vector.x=x/dValue;
		vector.y=y/dValue;
		vector.z=z/dValue;
		return vector;
	}
/*
		Vector3& Vector3::operator /= (float dValue){
			x/=dValue;y/=dValue;z/=dValue;
			return *this;
		}*/
	
	Matrix3 Vector3::operator * (Vector3& vector1){
		Matrix3 matrix3;
		matrix3.element[0][0]=x*vector1.x;
 		matrix3.element[0][1]=x*vector1.y;
 		matrix3.element[0][2]=x*vector1.z;
 		matrix3.element[1][0]=y*vector1.x;
 		matrix3.element[1][1]=y*vector1.y;
 		matrix3.element[1][2]=y*vector1.z;
 		matrix3.element[2][0]=z*vector1.x;
 		matrix3.element[2][1]=z*vector1.y;
 		matrix3.element[2][2]=z*vector1.z;
		return matrix3;
	}
	Vector3 Vector3::operator * (float dValue){
		Vector3	vector = *this;
		vector.x=x*dValue;
		vector.y=y*dValue;
		vector.z=z*dValue;
		return	vector;
	}
/*
		Vector3& Vector3::operator *= (float dValue){
			x*=dValue;
			y*=dValue;
			z*=dValue;
			return	*this;
		}*/
	
	Vector3& Vector3::operator += (Vector3& vector1){
		x+=vector1.x;
		y+=vector1.y;
		z+=vector1.z;
		return	*this;
	}
	float Vector3::amplitude(Vector3 vector){
		float result=(float)sqrt((vector.x)*(vector.x)+(vector.y)*(vector.y)+(vector.z)*(vector.z));
		return result;
	}
	float Vector3::Distance(Vector3 vector1,Vector3 vector2){return amplitude(vector1-vector2);}
	Vector3 Vector3::Direction(Vector3 vector1,Vector3 vector2){return (vector1-vector2)/Distance(vector1,vector2);}
	Vector3 Vector3::Normalize(Vector3 vector1){return vector1/vector.amplitude(vector1);}
	float Vector3::Max(Vector3 &vector){
		return MAX(MAX(vector.x,vector.y),vector.z);
	}
	Vector3 Vector3::mean(Vector3 *vertex,long nverts){
		Vector3 sum;
		for (int i=0;i<nverts;i++){
			sum+=vertex[i];
		}
		return sum/(float)nverts;
	};
	Vector3 Vector3::CrossOver(Vector3 vector1,Vector3 vector2){
		Vector3 result;
		result.x=vector1.y*vector2.z-vector2.y*vector1.z;
		result.y=vector1.x*vector2.z-vector2.x*vector1.z;
		result.z=vector1.x*vector2.y-vector2.x*vector1.y;
		return result;
	}

	Vector4 Vector3::Mat4MulV4(Matrix4 &m4,Vector4 &v4){
		Vector4 rv4;
		rv4.x=m4.element[0][0]*v4.x+m4.element[0][1]*v4.y+m4.element[0][2]*v4.z+m4.element[0][3]*v4.w;
		rv4.y=m4.element[1][0]*v4.x+m4.element[1][1]*v4.y+m4.element[1][2]*v4.z+m4.element[1][3]*v4.w;
		rv4.z=m4.element[2][0]*v4.x+m4.element[2][1]*v4.y+m4.element[2][2]*v4.z+m4.element[2][3]*v4.w;
		rv4.w=m4.element[3][0]*v4.x+m4.element[3][1]*v4.y+m4.element[3][2]*v4.z+m4.element[3][3]*v4.w;
		return rv4;
	}
	Vector3 Vector3::Mat3MulV3(Matrix3 &m3,Vector3 &v3){
		Vector3 rv3;
		rv3.x=m3.element[0][0]*v3.x+m3.element[0][1]*v3.y+m3.element[0][2]*v3.z;
		rv3.y=m3.element[1][0]*v3.x+m3.element[1][1]*v3.y+m3.element[1][2]*v3.z;
		rv3.z=m3.element[2][0]*v3.x+m3.element[2][1]*v3.y+m3.element[2][2]*v3.z;
		return rv3;
	}
	Vector3 Vector3::Rotate(Vector3 v,float a,float b,float r,float dx,float dy,float dz){
		Vector4 inv4(v);Matrix4 RotMat(a,b,r,dx,dy,dz);
		Vector4 outv4=vector.Mat4MulV4(RotMat,inv4);
		Vector3 outv3(outv4);
		return outv3;
	}


	
Vector4::Vector4(){x=y=z=0;w=1;}
Vector4::Vector4(Vector3 &p){x=p.x;y=p.y;z=p.z;w=1;}