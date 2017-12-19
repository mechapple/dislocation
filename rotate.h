void rotate(double R[3][3],double vec[3])
{
	double vec2[3];
	vec2[0] = R[0][0]*vec[0] + R[0][1]*vec[1] + R[0][2]*vec[2];
	vec2[1] = R[1][0]*vec[0] + R[1][1]*vec[1] + R[1][2]*vec[2];
	vec2[2] = R[2][0]*vec[0] + R[2][1]*vec[1] + R[2][2]*vec[2];
	
	vec[0] = vec2[0]; vec[1] = vec2[1]; vec[2] = vec2[2];
}

double angle_dot(double x1,double x2,double x3,double y1,double y2,double y3)
{
	double xmag = sqrt(x1*x1+x2*x2+x3*x3);
	double ymag = sqrt(y1*y1+y2*y2+y3*y3);
	
	double dot = x1*y1+x2*y2+x3*y3;
	double theta = acos(dot/(xmag*ymag));
	
	return theta;
}

void unit_cross(double x1,double x2,double x3,double y1,double y2,double y3,double n[3])
{
	n[0] = (x2*y3-x3*y2);
	n[1] = (x3*y1-x1*y3);
	n[2] = (x1*y2-x2*y1);
	
	double mag = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0] = n[0]/mag; n[1] = n[1]/mag; n[2] = n[2]/mag;
}
