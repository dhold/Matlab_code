function tmp = somefn(in1,in2,in3,in4)
%SOMEFN
%    TMP = SOMEFN(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 5.11.
%    15-Oct-2014 14:52:39

Fge1_1 = in1(1,:);
Fge1_2 = in1(5,:);
Fge1_3 = in1(9,:);
Fge1_4 = in1(13,:);
Fge1_5 = in1(17,:);
Fge1_6 = in1(21,:);
Fge1_7 = in1(25,:);
Fge1_8 = in1(29,:);
Gge1_1 = in2(1,:);
Gge1_2 = in2(5,:);
Gge1_3 = in2(9,:);
Gge1_4 = in2(13,:);
Gge1_5 = in2(17,:);
Gge1_6 = in2(21,:);
Gge1_7 = in2(25,:);
Gge1_8 = in2(29,:);
Gge2_1 = in2(2,:);
Gge2_2 = in2(6,:);
Gge2_3 = in2(10,:);
Gge2_4 = in2(14,:);
Gge2_5 = in2(18,:);
Gge2_6 = in2(22,:);
Gge2_7 = in2(26,:);
Gge2_8 = in2(30,:);
Gge3_1 = in2(3,:);
Gge3_2 = in2(7,:);
Gge3_3 = in2(11,:);
Gge3_4 = in2(15,:);
Gge3_5 = in2(19,:);
Gge3_6 = in2(23,:);
Gge3_7 = in2(27,:);
Gge3_8 = in2(31,:);
Gge4_1 = in2(4,:);
Gge4_2 = in2(8,:);
Gge4_3 = in2(12,:);
Gge4_4 = in2(16,:);
Gge4_5 = in2(20,:);
Gge4_6 = in2(24,:);
Gge4_7 = in2(28,:);
Gge4_8 = in2(32,:);
Hge1_1 = in3(1,:);
Hge1_2 = in3(5,:);
Hge1_3 = in3(9,:);
Hge1_4 = in3(13,:);
Hge1_5 = in3(17,:);
Hge1_6 = in3(21,:);
Hge1_7 = in3(25,:);
Hge1_8 = in3(29,:);
Hge2_1 = in3(2,:);
Hge2_2 = in3(6,:);
Hge2_3 = in3(10,:);
Hge2_4 = in3(14,:);
Hge2_5 = in3(18,:);
Hge2_6 = in3(22,:);
Hge2_7 = in3(26,:);
Hge2_8 = in3(30,:);
Hge3_1 = in3(3,:);
Hge3_2 = in3(7,:);
Hge3_3 = in3(11,:);
Hge3_4 = in3(15,:);
Hge3_5 = in3(19,:);
Hge3_6 = in3(23,:);
Hge3_7 = in3(27,:);
Hge3_8 = in3(31,:);
Hge4_1 = in3(4,:);
Hge4_2 = in3(8,:);
Hge4_3 = in3(12,:);
Hge4_4 = in3(16,:);
Hge4_5 = in3(20,:);
Hge4_6 = in3(24,:);
Hge4_7 = in3(28,:);
Hge4_8 = in3(32,:);
Ige1_1 = in4(1,:);
Ige1_2 = in4(5,:);
Ige1_3 = in4(9,:);
Ige1_4 = in4(13,:);
Ige1_5 = in4(17,:);
Ige1_6 = in4(21,:);
Ige1_7 = in4(25,:);
Ige1_8 = in4(29,:);
Ige2_1 = in4(2,:);
Ige2_2 = in4(6,:);
Ige2_3 = in4(10,:);
Ige2_4 = in4(14,:);
Ige2_5 = in4(18,:);
Ige2_6 = in4(22,:);
Ige2_7 = in4(26,:);
Ige2_8 = in4(30,:);
Ige3_1 = in4(3,:);
Ige3_2 = in4(7,:);
Ige3_3 = in4(11,:);
Ige3_4 = in4(15,:);
Ige3_5 = in4(19,:);
Ige3_6 = in4(23,:);
Ige3_7 = in4(27,:);
Ige3_8 = in4(31,:);
Ige4_1 = in4(4,:);
Ige4_2 = in4(8,:);
Ige4_3 = in4(12,:);
Ige4_4 = in4(16,:);
Ige4_5 = in4(20,:);
Ige4_6 = in4(24,:);
Ige4_7 = in4(28,:);
Ige4_8 = in4(32,:);
t2 = conj(Hge1_1);
t3 = Fge1_1.*t2;
t4 = conj(Hge1_2);
t5 = Fge1_2.*t4;
t6 = conj(Hge1_3);
t7 = Fge1_3.*t6;
t8 = conj(Hge1_4);
t9 = Fge1_4.*t8;
t10 = conj(Hge1_5);
t11 = Fge1_5.*t10;
t12 = conj(Hge1_6);
t13 = Fge1_6.*t12;
t14 = conj(Hge1_7);
t15 = Fge1_7.*t14;
t16 = conj(Hge1_8);
t17 = Fge1_8.*t16;
t18 = t3+t5+t7+t9+t11+t13+t15+t17;
t19 = conj(Hge2_1);
t20 = Fge1_1.*t19;
t21 = conj(Hge2_2);
t22 = Fge1_2.*t21;
t23 = conj(Hge2_3);
t24 = Fge1_3.*t23;
t25 = conj(Hge2_4);
t26 = Fge1_4.*t25;
t27 = conj(Hge2_5);
t28 = Fge1_5.*t27;
t29 = conj(Hge2_6);
t30 = Fge1_6.*t29;
t31 = conj(Hge2_7);
t32 = Fge1_7.*t31;
t33 = conj(Hge2_8);
t34 = Fge1_8.*t33;
t35 = t20+t22+t24+t26+t28+t30+t32+t34;
t36 = conj(Hge3_1);
t37 = Fge1_1.*t36;
t38 = conj(Hge3_2);
t39 = Fge1_2.*t38;
t40 = conj(Hge3_3);
t41 = Fge1_3.*t40;
t42 = conj(Hge3_4);
t43 = Fge1_4.*t42;
t44 = conj(Hge3_5);
t45 = Fge1_5.*t44;
t46 = conj(Hge3_6);
t47 = Fge1_6.*t46;
t48 = conj(Hge3_7);
t49 = Fge1_7.*t48;
t50 = conj(Hge3_8);
t51 = Fge1_8.*t50;
t52 = t37+t39+t41+t43+t45+t47+t49+t51;
t53 = conj(Hge4_1);
t54 = Fge1_1.*t53;
t55 = conj(Hge4_2);
t56 = Fge1_2.*t55;
t57 = conj(Hge4_3);
t58 = Fge1_3.*t57;
t59 = conj(Hge4_4);
t60 = Fge1_4.*t59;
t61 = conj(Hge4_5);
t62 = Fge1_5.*t61;
t63 = conj(Hge4_6);
t64 = Fge1_6.*t63;
t65 = conj(Hge4_7);
t66 = Fge1_7.*t65;
t67 = conj(Hge4_8);
t68 = Fge1_8.*t67;
t69 = t54+t56+t58+t60+t62+t64+t66+t68;
t70 = Gge1_1.*t2;
t71 = Gge1_2.*t4;
t72 = Gge1_3.*t6;
t73 = Gge1_4.*t8;
t74 = Gge1_5.*t10;
t75 = Gge1_6.*t12;
t76 = Gge1_7.*t14;
t77 = Gge1_8.*t16;
t78 = t70+t71+t72+t73+t74+t75+t76+t77;
t79 = Gge1_1.*t19;
t80 = Gge1_2.*t21;
t81 = Gge1_3.*t23;
t82 = Gge1_4.*t25;
t83 = Gge1_5.*t27;
t84 = Gge1_6.*t29;
t85 = Gge1_7.*t31;
t86 = Gge1_8.*t33;
t87 = t79+t80+t81+t82+t83+t84+t85+t86;
t88 = Gge1_1.*t36;
t89 = Gge1_2.*t38;
t90 = Gge1_3.*t40;
t91 = Gge1_4.*t42;
t92 = Gge1_5.*t44;
t93 = Gge1_6.*t46;
t94 = Gge1_7.*t48;
t95 = Gge1_8.*t50;
t96 = t88+t89+t90+t91+t92+t93+t94+t95;
t97 = Gge1_1.*t53;
t98 = Gge1_2.*t55;
t99 = Gge1_3.*t57;
t100 = Gge1_4.*t59;
t101 = Gge1_5.*t61;
t102 = Gge1_6.*t63;
t103 = Gge1_7.*t65;
t104 = Gge1_8.*t67;
t105 = t97+t98+t99+t100+t101+t102+t103+t104;
t106 = conj(Fge1_1);
t107 = conj(Fge1_2);
t108 = Ige1_1.*t2;
t109 = Ige1_2.*t4;
t110 = Ige1_3.*t6;
t111 = Ige1_4.*t8;
t112 = Ige1_5.*t10;
t113 = Ige1_6.*t12;
t114 = Ige1_7.*t14;
t115 = Ige1_8.*t16;
t116 = t108+t109+t110+t111+t112+t113+t114+t115;
t117 = Ige1_1.*t19;
t118 = Ige1_2.*t21;
t119 = Ige1_3.*t23;
t120 = Ige1_4.*t25;
t121 = Ige1_5.*t27;
t122 = Ige1_6.*t29;
t123 = Ige1_7.*t31;
t124 = Ige1_8.*t33;
t125 = t117+t118+t119+t120+t121+t122+t123+t124;
t126 = Ige1_1.*t36;
t127 = Ige1_2.*t38;
t128 = Ige1_3.*t40;
t129 = Ige1_4.*t42;
t130 = Ige1_5.*t44;
t131 = Ige1_6.*t46;
t132 = Ige1_7.*t48;
t133 = Ige1_8.*t50;
t134 = t126+t127+t128+t129+t130+t131+t132+t133;
t135 = Ige1_1.*t53;
t136 = Ige1_2.*t55;
t137 = Ige1_3.*t57;
t138 = Ige1_4.*t59;
t139 = Ige1_5.*t61;
t140 = Ige1_6.*t63;
t141 = Ige1_7.*t65;
t142 = Ige1_8.*t67;
t143 = t135+t136+t137+t138+t139+t140+t141+t142;
t144 = conj(Fge1_3);
t145 = conj(Fge1_4);
t146 = conj(Fge1_5);
t147 = conj(Fge1_6);
t148 = conj(Fge1_7);
t149 = conj(Fge1_8);
t150 = conj(Gge1_1);
t151 = conj(Gge1_2);
t152 = conj(Gge1_3);
t153 = conj(Gge1_4);
t154 = conj(Gge1_5);
t155 = conj(Gge1_6);
t156 = conj(Gge1_7);
t157 = conj(Gge1_8);
t158 = Fge1_1.*t150;
t159 = Fge1_2.*t151;
t160 = Fge1_3.*t152;
t161 = Fge1_4.*t153;
t162 = Fge1_5.*t154;
t163 = Fge1_6.*t155;
t164 = Fge1_7.*t156;
t165 = Fge1_8.*t157;
t166 = t158+t159+t160+t161+t162+t163+t164+t165;
t167 = conj(Gge2_1);
t168 = Fge1_1.*t167;
t169 = conj(Gge2_2);
t170 = Fge1_2.*t169;
t171 = conj(Gge2_3);
t172 = Fge1_3.*t171;
t173 = conj(Gge2_4);
t174 = Fge1_4.*t173;
t175 = conj(Gge2_5);
t176 = Fge1_5.*t175;
t177 = conj(Gge2_6);
t178 = Fge1_6.*t177;
t179 = conj(Gge2_7);
t180 = Fge1_7.*t179;
t181 = conj(Gge2_8);
t182 = Fge1_8.*t181;
t183 = t168+t170+t172+t174+t176+t178+t180+t182;
t184 = conj(Gge3_1);
t185 = Fge1_1.*t184;
t186 = conj(Gge3_2);
t187 = Fge1_2.*t186;
t188 = conj(Gge3_3);
t189 = Fge1_3.*t188;
t190 = conj(Gge3_4);
t191 = Fge1_4.*t190;
t192 = conj(Gge3_5);
t193 = Fge1_5.*t192;
t194 = conj(Gge3_6);
t195 = Fge1_6.*t194;
t196 = conj(Gge3_7);
t197 = Fge1_7.*t196;
t198 = conj(Gge3_8);
t199 = Fge1_8.*t198;
t200 = t185+t187+t189+t191+t193+t195+t197+t199;
t201 = conj(Gge4_1);
t202 = Fge1_1.*t201;
t203 = conj(Gge4_2);
t204 = Fge1_2.*t203;
t205 = conj(Gge4_3);
t206 = Fge1_3.*t205;
t207 = conj(Gge4_4);
t208 = Fge1_4.*t207;
t209 = conj(Gge4_5);
t210 = Fge1_5.*t209;
t211 = conj(Gge4_6);
t212 = Fge1_6.*t211;
t213 = conj(Gge4_7);
t214 = Fge1_7.*t213;
t215 = conj(Gge4_8);
t216 = Fge1_8.*t215;
t217 = t202+t204+t206+t208+t210+t212+t214+t216;
tmp = imag(t106.*(Gge1_1.*t116+Gge2_1.*t125+Gge3_1.*t134+Gge4_1.*t143)).*2.0+imag(t107.*(Gge1_2.*t116+Gge2_2.*t125+Gge3_2.*t134+Gge4_2.*t143)).*2.0+imag(t144.*(Gge1_3.*t116+Gge2_3.*t125+Gge3_3.*t134+Gge4_3.*t143)).*2.0+imag(t145.*(Gge1_4.*t116+Gge2_4.*t125+Gge3_4.*t134+Gge4_4.*t143)).*2.0+imag(t146.*(Gge1_5.*t116+Gge2_5.*t125+Gge3_5.*t134+Gge4_5.*t143)).*2.0+imag(t147.*(Gge1_6.*t116+Gge2_6.*t125+Gge3_6.*t134+Gge4_6.*t143)).*2.0+imag(t148.*(Gge1_7.*t116+Gge2_7.*t125+Gge3_7.*t134+Gge4_7.*t143)).*2.0+imag(t149.*(Gge1_8.*t116+Gge2_8.*t125+Gge3_8.*t134+Gge4_8.*t143)).*2.0+imag(t150.*(Ige1_1.*t18+Ige2_1.*t35+Ige3_1.*t52+Ige4_1.*t69)).*2.0+imag(t151.*(Ige1_2.*t18+Ige2_2.*t35+Ige3_2.*t52+Ige4_2.*t69)).*2.0+imag(t152.*(Ige1_3.*t18+Ige2_3.*t35+Ige3_3.*t52+Ige4_3.*t69)).*2.0+imag(t153.*(Ige1_4.*t18+Ige2_4.*t35+Ige3_4.*t52+Ige4_4.*t69)).*2.0+imag(t154.*(Ige1_5.*t18+Ige2_5.*t35+Ige3_5.*t52+Ige4_5.*t69)).*2.0+imag(t155.*(Ige1_6.*t18+Ige2_6.*t35+Ige3_6.*t52+Ige4_6.*t69)).*2.0+imag(t156.*(Ige1_7.*t18+Ige2_7.*t35+Ige3_7.*t52+Ige4_7.*t69)).*2.0+imag(t157.*(Ige1_8.*t18+Ige2_8.*t35+Ige3_8.*t52+Ige4_8.*t69)).*2.0+imag(t106.*(Ige1_1.*t78+Ige2_1.*t87+Ige3_1.*t96+Ige4_1.*t105)).*2.0+imag(t107.*(Ige1_2.*t78+Ige2_2.*t87+Ige3_2.*t96+Ige4_2.*t105)).*2.0+imag(t144.*(Ige1_3.*t78+Ige2_3.*t87+Ige3_3.*t96+Ige4_3.*t105)).*2.0+imag(t145.*(Ige1_4.*t78+Ige2_4.*t87+Ige3_4.*t96+Ige4_4.*t105)).*2.0+imag(t146.*(Ige1_5.*t78+Ige2_5.*t87+Ige3_5.*t96+Ige4_5.*t105)).*2.0+imag(t147.*(Ige1_6.*t78+Ige2_6.*t87+Ige3_6.*t96+Ige4_6.*t105)).*2.0+imag(t148.*(Ige1_7.*t78+Ige2_7.*t87+Ige3_7.*t96+Ige4_7.*t105)).*2.0+imag(t149.*(Ige1_8.*t78+Ige2_8.*t87+Ige3_8.*t96+Ige4_8.*t105)).*2.0+imag(t2.*(Ige1_1.*t166+Ige2_1.*t183+Ige3_1.*t200+Ige4_1.*t217)).*2.0+imag(t4.*(Ige1_2.*t166+Ige2_2.*t183+Ige3_2.*t200+Ige4_2.*t217)).*2.0+imag(t6.*(Ige1_3.*t166+Ige2_3.*t183+Ige3_3.*t200+Ige4_3.*t217)).*2.0+imag(t8.*(Ige1_4.*t166+Ige2_4.*t183+Ige3_4.*t200+Ige4_4.*t217)).*2.0+imag(t10.*(Ige1_5.*t166+Ige2_5.*t183+Ige3_5.*t200+Ige4_5.*t217)).*2.0+imag(t12.*(Ige1_6.*t166+Ige2_6.*t183+Ige3_6.*t200+Ige4_6.*t217)).*2.0+imag(t14.*(Ige1_7.*t166+Ige2_7.*t183+Ige3_7.*t200+Ige4_7.*t217)).*2.0+imag(t16.*(Ige1_8.*t166+Ige2_8.*t183+Ige3_8.*t200+Ige4_8.*t217)).*2.0;
