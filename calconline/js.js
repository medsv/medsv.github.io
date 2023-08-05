function SaturationPressure(t){
var theta, A, b , C;
//var SaturationPressure;
    var n =[1167.0521452767, -724213.16703206, -17.073846940092, 12020.82470247, -3232555.0322333, 14.91510861353, -4823.2657361591, 405113.40542057, -0.23855557567849, 650.17534844798];
    t = t + 273.15;
    theta = t + n[8] / (t - n[9]);
    A = theta * theta + n[0] * theta + n[1];
    b = n[2] * theta * theta + n[3] * theta + n[4];
    C = n[5] * theta * theta + n[6] * theta + n[7];
//    SaturationPressure = (2 * C / (-b + (b * b - 4 * A * C) ** 0.5)) ** 4 * 1000000;
	return (2 * C / (-b + (b * b - 4 * A * C) ** 0.5)) ** 4 * 1000000;
//document.write(SaturationPressure);
}

function  SaturationTemperature(P){
    var beta, E, F, G, D;
    var n = [1167.0521452767, -724213.16703206, -17.073846940092, 12020.82470247, -3232555.0322333, 14.91510861353, -4823.2657361591, 405113.40542057, -0.23855557567849, 650.17534844798];
    beta = (P / 1000000) ** 0.25;
    E = beta * beta + n[2] * beta + n[5];
    F = n[0] * beta * beta + n[3] * beta + n[6];
    G = n[1] * beta * beta + n[4] * beta + n[7];
    D = 2 * G / (-F - (F * F - 4 * E * G) ** 0.5);
    return (n[9] + D - ((n[9] + D) ** 2 - 4 * (n[8] + n[9] * D)) ** 0.5) / 2 - 273.15;
}

function enthalpy(t,d){
	return 1.006*t+d*(1.84*t+2501)/ 1000;
}

function humidity (p, t, RH){
	var p_s, p_st, p_da, ro_v, ro_wd, d;
    var mu_st = 18.;  // Молярная масса водяного пара, г/моль
    var mu_da = 29.;  // Молярная масса сухого воздуха, г/моль
    var R = 8.314  // Универсальная газовая постоянная, Дж/(моль*К) 
	p_s=SaturationPressure(t);
	p_st=RH*p_s; 
//	ro_v=P_v/461.495/(t+273.15);
    p_da=p-p_st;
//	ro_da=P_da/287.058/(t+273.15);
//	d=ro_v/ro_da;
//	var ro_wd=ro_v+ro_da;
    d = 1000 * p_st / (p - p_st) * mu_st / mu_da;
    ro_wd = (mu_st * p_st + mu_da * p_da) / R / (t + 273.15) / 1000.;
    
	return [d, ro_wd];
}

function wetTermometer (P,t_d, RH){
//_v - vapor
//_da - dry air
//_wa - wet air
//_d - dry
//_w - wet
//_dp - dew point
//_s - saturated
//ro - density
	var t_w, t_dp, Psat, P_v, d;
	var eps=0.001;
	var answer =[];
	Psat=SaturationPressure(t_d);
	P_v=RH*Psat; 
	t_dp=SaturationTemperature(P_v);
	answer=humidity(P, t_d, RH);
	d=answer[0];
	
	var x0, x1, y0, y1, yy, dydx;
	yy=enthalpy(t_d, d); //цель
	x0=0.75*t_d;	//первое приближение
	answer=humidity(P, x0, 1.);
	d=answer[0];
	y0=enthalpy(x0, d);
	x1=x0+1;
	answer=humidity(P, x1, 1.);
	d=answer[0];
	y1=enthalpy(x1, d);
	dydx=y1-y0;
	
	do {
	y1=yy;
	x1=x0+(y1-y0)/dydx;
	answer=humidity(P,x1,1.0);
	d=answer[0];
	
	y1=enthalpy(x1, d);
	dydx=(y1-y0)/(x1-x0);
	x0=x1;
	y0=y1;
	} while (Math.abs(y1-yy)>eps);

	return [x1, t_dp];

}
