#ifndef VECTOR_H_FILE
#define VECTOR_H_FILE

//作者:  黄海波，张嘉骞，中国科学技术大学

template<class myType>
class vector3 {

public:
	vector3(myType x, myType y, myType z);
	vector3();
	myType dot(vector3<myType> other);
	vector3<myType> cross(vector3<myType> other);
	myType norm_square();
	myType x;
	myType y;
	myType z;
	myType get_x()	const
	{
		return x;
	}
	myType get_y()	const
	{
		return y;
	}
	myType get_z()	const
	{
		return z;
	}
	vector3<myType> operator- (const vector3<myType>& rhs) const {
		return vector3<myType>(x - rhs.x, y - rhs.y, z - rhs.z);
	}
	vector3<myType> operator+ (const vector3<myType>& rhs) const {
		return vector3<myType>(x + rhs.x, y + rhs.y, z + rhs.z);
	}
	vector3<myType> operator/ (const myType& rhs) const {
		return vector3<myType>(x/rhs, y/rhs, z/rhs);
	}
};

#endif
