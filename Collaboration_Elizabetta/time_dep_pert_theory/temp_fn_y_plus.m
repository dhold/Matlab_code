function out =temp_fn_y_plus(om,om1,om2,om_u,t_sep) 
out =1.0.*exp(t_sep.*(om1 + om2).*(0.0 - 1.0.*i)).*(((om1 - 1.0.*om + om2).*(0.0000000005537407725647013.*i + 0.0000002276659150255968))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(76.16059769692698.*i + 16707.75673560291) + om2.*(45.69635861815619.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1392.104793972758 - 508989.0956633501.*i))) + (om1.*(0.00000000009275788318708888.*i + 0.000000002970526973009925))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(17332.06223846845 - 122.9980484076083.*i) + om2.*(17332.06223846845 - 73.79882904456498.*i) + om1.*om2 + om1.^2 - (852723.9320843287.*i + 3630.844778899285))) + (om2.*(- 0.0000000004202263172708945.*i - 0.000000001249636162341404))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(17332.06223846845 - 122.9980484076083.*i) + om2.*(17332.06223846845 - 73.79882904456498.*i) + om1.*om2 + om1.^2 - (852723.9320843287.*i + 3630.844778899285))) + ((om1 - 1.0.*om + om2).*(0.000000000182252721769408.*i + 0.00000002790393366564496))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(16707.75673560291 - 113.630558265472.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (822008.5886992921.*i + 3169.971576515177))) + (om1.*(- 0.0000000007163758346636341.*i - 0.0000001978855184684595))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(85.52808783906324.*i + 17332.06223846845) + om2.*(55.06384876029245.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1677.478253230826 - 528008.087760838.*i))) + (om2.*(0.0000000002716461518966096.*i + 0.00000002904488354081542))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(85.52808783906324.*i + 17332.06223846845) + om2.*(55.06384876029245.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1677.478253230826 - 528008.087760838.*i))) + ((om1 - 1.0.*om + om2).*(0.0000000003134102032345504.*i + 0.0000000007040846071588818))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(15459.14572987184 - 94.89557798119951.*i) + om2.*(16083.45123273737 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (606255.0106500527.*i + 10043180.39798112))) + (om1.*(0.0000000007163758346636341.*i + 0.0000001978855184684595))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(- 94.89557798119951.*i - 17956.36774133398) + om2.*(- 55.06384876029245.*i - 17332.06223846845) + om1.*om2 + om1.^2 + (724742.6737125896.*i + 10818308.54317013))) + (om2.*(- 0.0000000002716461518966096.*i - 0.00000002904488354081542))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(- 94.89557798119951.*i - 17956.36774133398) + om2.*(- 55.06384876029245.*i - 17332.06223846845) + om1.*om2 + om1.^2 + (724742.6737125896.*i + 10818308.54317013))) + ((om1 - 1.0.*om + om2).*(0.00000002908519946984528 - 0.0000000001143464823833147.*i))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(113.630558265472.*i + 16707.75673560291) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3169.971576515177 - 822008.5886992921.*i))) + (om1.*(0.0000000004736445170227625.*i - 0.000000001442104960926611))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(18580.67324419952 - 104.2630681233358.*i) + om2.*(17956.36774133398 - 64.43133890242872.*i) + om1.*om2 + om1.^2 + (11207692.7807475 - 755458.0170976262.*i))) + (om2.*(0.0000000003707402041576675 - 0.0000000001078491224439891.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(18580.67324419952 - 104.2630681233358.*i) + om2.*(17956.36774133398 - 64.43133890242872.*i) + om1.*om2 + om1.^2 + (11207692.7807475 - 755458.0170976262.*i))) + (om1.*(0.0000000005954551146454197.*i + 0.0000003942021366573875))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(76.16059769692698.*i + 16707.75673560291) + om2.*(45.69635861815619.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1392.104793972758 - 508989.0956633501.*i))) + (om2.*(- 0.0000000003585456101852984.*i - 0.0000002276713260370336))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(76.16059769692698.*i + 16707.75673560291) + om2.*(45.69635861815619.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1392.104793972758 - 508989.0956633501.*i))) + ((om1 - 1.0.*om + om2).*(0.0000000004998984867928307.*i - 0.0000002277318339358096))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(16707.75673560291 - 76.16059769692698.*i) + om2.*(16707.75673560291 - 45.69635861815619.*i) + om1.*om2 + om1.^2 - (508989.0956633501.*i + 1392.104793972758))) + (om1.*(- 0.0000000003178001393999088.*i - 0.000000002095018642718141))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(16707.75673560291 - 113.630558265472.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (822008.5886992921.*i + 3169.971576515177))) + (om2.*(0.000000000532971342409152.*i + 0.00000002871795531103577))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(16707.75673560291 - 113.630558265472.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (822008.5886992921.*i + 3169.971576515177))) + (om1.*(0.0000000001948172358816735.*i - 0.000000195721865291089))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(113.630558265472.*i + 16707.75673560291) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3169.971576515177 - 822008.5886992921.*i))) + (om2.*(0.00000002803981290599838 - 0.00000000006086649043893078.*i))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(113.630558265472.*i + 16707.75673560291) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3169.971576515177 - 822008.5886992921.*i))) + ((om1 - 1.0.*om + om2).*(- 0.0000000005365002362780465.*i - 0.00000002759894143120577))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(17332.06223846845 - 85.52808783906324.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (528008.087760838.*i + 1677.478253230826))) + (om1.*(0.0000000001948172358816735.*i + 0.000000195721865291089))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(16083.45123273737 - 104.2630681233358.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (625274.0027475406.*i + 10433310.88222013))) + (om2.*(- 0.00000000006086649043893078.*i - 0.00000002803981290599838))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(16083.45123273737 - 104.2630681233358.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (625274.0027475406.*i + 10433310.88222013))) + (om1.*(- 0.0000000004736445170227625.*i - 0.000000001442104960926611))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(- 104.2630681233358.*i - 18580.67324419952) + om2.*(- 64.43133890242872.*i - 17956.36774133398) + om1.*om2 + om1.^2 + (755458.0170976262.*i + 11207692.7807475))) + (om2.*(0.0000000001078491224439891.*i + 0.0000000003707402041576675))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(- 104.2630681233358.*i - 18580.67324419952) + om2.*(- 64.43133890242872.*i - 17956.36774133398) + om1.*om2 + om1.^2 + (755458.0170976262.*i + 11207692.7807475))) + ((om1 - 1.0.*om + om2).*(0.00000000005355378882752603.*i - 0.000000001715034547502333))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(122.9980484076083.*i + 17332.06223846845) + om2.*(73.79882904456498.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3630.844778899285 - 852723.9320843287.*i))) + ((om1 - 1.0.*om + om2).*(0.0000000001114725148867249.*i - 0.0000000007804759174949449))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(94.89557798119951.*i + 15459.14572987184) + om2.*(55.06384876029245.*i + 16083.45123273737) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10043180.39798112 - 606255.0106500527.*i))) + (om1.*(0.0000000009591071523045057.*i - 0.0000003943289319759923))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(16707.75673560291 - 76.16059769692698.*i) + om2.*(16707.75673560291 - 45.69635861815619.*i) + om1.*om2 + om1.^2 - (508989.0956633501.*i + 1392.104793972758))) + (om2.*(0.0000002274003639950666 - 0.0000000003101178656670951.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(16707.75673560291 - 76.16059769692698.*i) + om2.*(16707.75673560291 - 45.69635861815619.*i) + om1.*om2 + om1.^2 - (508989.0956633501.*i + 1392.104793972758))) + (om1.*(0.0000000005428423956127287.*i - 0.000000001219510312426357))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(94.89557798119951.*i + 15459.14572987184) + om2.*(55.06384876029245.*i + 16083.45123273737) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10043180.39798112 - 606255.0106500527.*i))) + (om2.*(0.0000000004508440612794897 - 0.0000000007710416832295444.*i))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(94.89557798119951.*i + 15459.14572987184) + om2.*(55.06384876029245.*i + 16083.45123273737) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10043180.39798112 - 606255.0106500527.*i))) + ((om1 - 1.0.*om + om2).*(0.00000002941651349098637 - 0.000000000482735294631715.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(17956.36774133398 - 94.89557798119951.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 + (10818308.54317013 - 724742.6737125896.*i))) + (om1.*(0.0000000004634705225378483.*i + 0.000000001110853722912134))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(15459.14572987184 - 94.89557798119951.*i) + om2.*(16083.45123273737 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (606255.0106500527.*i + 10043180.39798112))) + (om2.*(- 0.0000000002110539904719353.*i - 0.0000000007427509034929421))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(15459.14572987184 - 94.89557798119951.*i) + om2.*(16083.45123273737 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (606255.0106500527.*i + 10043180.39798112))) + ((om1 - 1.0.*om + om2).*(0.00000002908519946984528 - 0.0000000001143464823833147.*i))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(104.2630681233358.*i + 16083.45123273737) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10433310.88222013 - 625274.0027475406.*i))) + (om1.*(0.0000001978855184684595 - 0.0000000007163758346636341.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(17332.06223846845 - 85.52808783906324.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (528008.087760838.*i + 1677.478253230826))) + (om2.*(0.0000000002716461518966096.*i - 0.00000002904488354081542))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(17332.06223846845 - 85.52808783906324.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (528008.087760838.*i + 1677.478253230826))) + (om1.*(0.00000000008244215340618773.*i - 0.0000000020832593705591))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(122.9980484076083.*i + 17332.06223846845) + om2.*(73.79882904456498.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3630.844778899285 - 852723.9320843287.*i))) + (om2.*(0.000000001340551445466435 - 0.00000000008202611053651234.*i))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(122.9980484076083.*i + 17332.06223846845) + om2.*(73.79882904456498.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3630.844778899285 - 852723.9320843287.*i))) + (om1.*(0.00000000002304386273861773 - 0.0000000008083919294656287.*i))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(17956.36774133398 - 94.89557798119951.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 + (10818308.54317013 - 724742.6737125896.*i))) + (om2.*(0.0000000003798154597319033.*i + 0.00000002776552824962252))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(17956.36774133398 - 94.89557798119951.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 + (10818308.54317013 - 724742.6737125896.*i))) + ((om1 - 1.0.*om + om2).*(0.00000000005355378882752603.*i + 0.000000001715034547502333))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(17332.06223846845 - 122.9980484076083.*i) + om2.*(17332.06223846845 - 73.79882904456498.*i) + om1.*om2 + om1.^2 - (852723.9320843287.*i + 3630.844778899285))) + ((om1 - 1.0.*om + om2).*(- 0.000000000482735294631715.*i - 0.00000002941651349098637))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(85.52808783906324.*i + 17332.06223846845) + om2.*(55.06384876029245.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1677.478253230826 - 528008.087760838.*i))) + ((om1 - 1.0.*om + om2).*(0.000000000482735294631715.*i + 0.00000002941651349098637))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(- 94.89557798119951.*i - 17956.36774133398) + om2.*(- 55.06384876029245.*i - 17332.06223846845) + om1.*om2 + om1.^2 + (724742.6737125896.*i + 10818308.54317013))) + (om1.*(0.0000000008083919294656287.*i + 0.00000000002304386273861773))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(- 94.89557798119951.*i - 17956.36774133398) + om2.*(- 55.06384876029245.*i - 17332.06223846845) + om1.*om2 + om1.^2 + (724742.6737125896.*i + 10818308.54317013))) + (om2.*(0.00000002776552824962252 - 0.0000000003798154597319033.*i))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(- 94.89557798119951.*i - 17956.36774133398) + om2.*(- 55.06384876029245.*i - 17332.06223846845) + om1.*om2 + om1.^2 + (724742.6737125896.*i + 10818308.54317013))) + ((om1 - 1.0.*om + om2).*(0.0000000002734587894032816.*i - 0.0000000008325996873906735))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(18580.67324419952 - 104.2630681233358.*i) + om2.*(17956.36774133398 - 64.43133890242872.*i) + om1.*om2 + om1.^2 + (11207692.7807475 - 755458.0170976262.*i))) + ((om1 - 1.0.*om + om2).*(0.0000000004998984867928307.*i + 0.0000002277318339358096))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(76.16059769692698.*i + 16707.75673560291) + om2.*(45.69635861815619.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1392.104793972758 - 508989.0956633501.*i))) + (om1.*(0.0000000001948172358816735.*i - 0.000000195721865291089))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(104.2630681233358.*i + 16083.45123273737) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10433310.88222013 - 625274.0027475406.*i))) + (om2.*(0.00000002803981290599838 - 0.00000000006086649043893078.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(104.2630681233358.*i + 16083.45123273737) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10433310.88222013 - 625274.0027475406.*i))) + (om1.*(0.00000000008578159657845464.*i - 0.0000000004817947437389991))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(18580.67324419952 - 104.2630681233358.*i) + om2.*(17956.36774133398 - 64.43133890242872.*i) + om1.*om2 + om1.^2 + (11207692.7807475 - 755458.0170976262.*i))) + (om2.*(0.00000000001372777260780094.*i + 0.0000000003773131407215033))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(18580.67324419952 - 104.2630681233358.*i) + om2.*(17956.36774133398 - 64.43133890242872.*i) + om1.*om2 + om1.^2 + (11207692.7807475 - 755458.0170976262.*i))) + ((om1 - 1.0.*om + om2).*(- 0.0000000001143464823833147.*i - 0.00000002908519946984528))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(16707.75673560291 - 113.630558265472.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (822008.5886992921.*i + 3169.971576515177))) + (om1.*(- 0.0000000008083919294656287.*i - 0.00000000002304386273861773))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(85.52808783906324.*i + 17332.06223846845) + om2.*(55.06384876029245.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1677.478253230826 - 528008.087760838.*i))) + (om2.*(0.0000000003798154597319033.*i - 0.00000002776552824962252))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(85.52808783906324.*i + 17332.06223846845) + om2.*(55.06384876029245.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1677.478253230826 - 528008.087760838.*i))) + ((om1 - 1.0.*om + om2).*(0.000000000182252721769408.*i - 0.00000002790393366564496))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(113.630558265472.*i + 16707.75673560291) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3169.971576515177 - 822008.5886992921.*i))) + (om1.*(0.00000000008244215340618773.*i + 0.0000000020832593705591))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(17332.06223846845 - 122.9980484076083.*i) + om2.*(17332.06223846845 - 73.79882904456498.*i) + om1.*om2 + om1.^2 - (852723.9320843287.*i + 3630.844778899285))) + (om2.*(- 0.00000000008202611053651234.*i - 0.000000001340551445466435))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(17332.06223846845 - 122.9980484076083.*i) + om2.*(17332.06223846845 - 73.79882904456498.*i) + om1.*om2 + om1.^2 - (852723.9320843287.*i + 3630.844778899285))) + ((om1 - 1.0.*om + om2).*(0.000000000182252721769408.*i + 0.00000002790393366564496))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(16083.45123273737 - 104.2630681233358.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (625274.0027475406.*i + 10433310.88222013))) + ((om1 - 1.0.*om + om2).*(- 0.0000000002734587894032816.*i - 0.0000000008325996873906735))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(- 104.2630681233358.*i - 18580.67324419952) + om2.*(- 64.43133890242872.*i - 17956.36774133398) + om1.*om2 + om1.^2 + (755458.0170976262.*i + 11207692.7807475))) + (om1.*(- 0.0000000003178001393999088.*i - 0.000000002095018642718141))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(16083.45123273737 - 104.2630681233358.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (625274.0027475406.*i + 10433310.88222013))) + (om2.*(0.000000000532971342409152.*i + 0.00000002871795531103577))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(16083.45123273737 - 104.2630681233358.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (625274.0027475406.*i + 10433310.88222013))) + (om1.*(- 0.00000000008578159657845464.*i - 0.0000000004817947437389991))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(- 104.2630681233358.*i - 18580.67324419952) + om2.*(- 64.43133890242872.*i - 17956.36774133398) + om1.*om2 + om1.^2 + (755458.0170976262.*i + 11207692.7807475))) + (om2.*(0.0000000003773131407215033 - 0.00000000001372777260780094.*i))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(- 104.2630681233358.*i - 18580.67324419952) + om2.*(- 64.43133890242872.*i - 17956.36774133398) + om1.*om2 + om1.^2 + (755458.0170976262.*i + 11207692.7807475))) + ((om1 - 1.0.*om + om2).*(0.0000000005537407725647013.*i - 0.0000002276659150255968))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(16707.75673560291 - 76.16059769692698.*i) + om2.*(16707.75673560291 - 45.69635861815619.*i) + om1.*om2 + om1.^2 - (508989.0956633501.*i + 1392.104793972758))) + (om1.*(0.0000000009591071523045057.*i + 0.0000003943289319759923))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(76.16059769692698.*i + 16707.75673560291) + om2.*(45.69635861815619.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1392.104793972758 - 508989.0956633501.*i))) + (om2.*(- 0.0000000003101178656670951.*i - 0.0000002274003639950666))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(76.16059769692698.*i + 16707.75673560291) + om2.*(45.69635861815619.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1392.104793972758 - 508989.0956633501.*i))) + ((om1 - 1.0.*om + om2).*(0.0000000003134102032345504.*i - 0.0000000007040846071588818))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(94.89557798119951.*i + 15459.14572987184) + om2.*(55.06384876029245.*i + 16083.45123273737) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10043180.39798112 - 606255.0106500527.*i))) + ((om1 - 1.0.*om + om2).*(0.0000000001114725148867249.*i + 0.0000000007804759174949449))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(15459.14572987184 - 94.89557798119951.*i) + om2.*(16083.45123273737 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (606255.0106500527.*i + 10043180.39798112))) + (om1.*(0.0000000001948172358816735.*i + 0.000000195721865291089))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(16707.75673560291 - 113.630558265472.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (822008.5886992921.*i + 3169.971576515177))) + (om2.*(- 0.00000000006086649043893078.*i - 0.00000002803981290599838))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(16707.75673560291 - 113.630558265472.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (822008.5886992921.*i + 3169.971576515177))) + ((om1 - 1.0.*om + om2).*(0.00000002941651349098637 - 0.000000000482735294631715.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(17332.06223846845 - 85.52808783906324.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (528008.087760838.*i + 1677.478253230826))) + ((om1 - 1.0.*om + om2).*(- 0.00000000006116487653923178.*i - 0.000000001574584912171271))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(122.9980484076083.*i + 17332.06223846845) + om2.*(73.79882904456498.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3630.844778899285 - 852723.9320843287.*i))) + (om1.*(0.0000000005428423956127287.*i + 0.000000001219510312426357))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(15459.14572987184 - 94.89557798119951.*i) + om2.*(16083.45123273737 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (606255.0106500527.*i + 10043180.39798112))) + (om2.*(- 0.0000000007710416832295444.*i - 0.0000000004508440612794897))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(15459.14572987184 - 94.89557798119951.*i) + om2.*(16083.45123273737 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (606255.0106500527.*i + 10043180.39798112))) + (om1.*(0.000000002095018642718141 - 0.0000000003178001393999088.*i))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(113.630558265472.*i + 16707.75673560291) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3169.971576515177 - 822008.5886992921.*i))) + (om2.*(0.000000000532971342409152.*i - 0.00000002871795531103577))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(113.630558265472.*i + 16707.75673560291) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3169.971576515177 - 822008.5886992921.*i))) + ((om1 - 1.0.*om + om2).*(- 0.0000000005365002362780465.*i - 0.00000002759894143120577))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(17956.36774133398 - 94.89557798119951.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 + (10818308.54317013 - 724742.6737125896.*i))) + (om1.*(0.0000000005954551146454197.*i - 0.0000003942021366573875))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(16707.75673560291 - 76.16059769692698.*i) + om2.*(16707.75673560291 - 45.69635861815619.*i) + om1.*om2 + om1.^2 - (508989.0956633501.*i + 1392.104793972758))) + (om2.*(0.0000002276713260370336 - 0.0000000003585456101852984.*i))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(16707.75673560291 - 76.16059769692698.*i) + om2.*(16707.75673560291 - 45.69635861815619.*i) + om1.*om2 + om1.^2 - (508989.0956633501.*i + 1392.104793972758))) + ((om1 - 1.0.*om + om2).*(0.0000000005365002362780465.*i - 0.00000002759894143120577))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(- 94.89557798119951.*i - 17956.36774133398) + om2.*(- 55.06384876029245.*i - 17332.06223846845) + om1.*om2 + om1.^2 + (724742.6737125896.*i + 10818308.54317013))) + (om1.*(0.00000000002304386273861773 - 0.0000000008083919294656287.*i))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(17332.06223846845 - 85.52808783906324.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (528008.087760838.*i + 1677.478253230826))) + (om2.*(0.0000000003798154597319033.*i + 0.00000002776552824962252))./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(17332.06223846845 - 85.52808783906324.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (528008.087760838.*i + 1677.478253230826))) + (om1.*(0.00000000009275788318708888.*i - 0.000000002970526973009925))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(122.9980484076083.*i + 17332.06223846845) + om2.*(73.79882904456498.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3630.844778899285 - 852723.9320843287.*i))) + (om2.*(0.000000001249636162341404 - 0.0000000004202263172708945.*i))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(122.9980484076083.*i + 17332.06223846845) + om2.*(73.79882904456498.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3630.844778899285 - 852723.9320843287.*i))) + (om1.*(0.0000000004634705225378483.*i - 0.000000001110853722912134))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(94.89557798119951.*i + 15459.14572987184) + om2.*(55.06384876029245.*i + 16083.45123273737) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10043180.39798112 - 606255.0106500527.*i))) + (om2.*(0.0000000007427509034929421 - 0.0000000002110539904719353.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(94.89557798119951.*i + 15459.14572987184) + om2.*(55.06384876029245.*i + 16083.45123273737) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10043180.39798112 - 606255.0106500527.*i))) + ((om1 - 1.0.*om + om2).*(0.000000000182252721769408.*i - 0.00000002790393366564496))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(104.2630681233358.*i + 16083.45123273737) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10433310.88222013 - 625274.0027475406.*i))) + ((om1 - 1.0.*om + om2).*(0.0000000001582889038769531.*i - 0.0000000006499788787790134))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(18580.67324419952 - 104.2630681233358.*i) + om2.*(17956.36774133398 - 64.43133890242872.*i) + om1.*om2 + om1.^2 + (11207692.7807475 - 755458.0170976262.*i))) + ((om1 - 1.0.*om + om2).*(0.00000002759894143120577 - 0.0000000005365002362780465.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(85.52808783906324.*i + 17332.06223846845) + om2.*(55.06384876029245.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1677.478253230826 - 528008.087760838.*i))) + ((om1 - 1.0.*om + om2).*(0.000000001574584912171271 - 0.00000000006116487653923178.*i))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(17332.06223846845 - 122.9980484076083.*i) + om2.*(17332.06223846845 - 73.79882904456498.*i) + om1.*om2 + om1.^2 - (852723.9320843287.*i + 3630.844778899285))) + (om1.*(0.0000001978855184684595 - 0.0000000007163758346636341.*i))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(17956.36774133398 - 94.89557798119951.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 + (10818308.54317013 - 724742.6737125896.*i))) + (om2.*(0.0000000002716461518966096.*i - 0.00000002904488354081542))./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(17956.36774133398 - 94.89557798119951.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 + (10818308.54317013 - 724742.6737125896.*i))) + (om1.*(0.000000002095018642718141 - 0.0000000003178001393999088.*i))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(104.2630681233358.*i + 16083.45123273737) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10433310.88222013 - 625274.0027475406.*i))) + (om2.*(0.000000000532971342409152.*i - 0.00000002871795531103577))./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(104.2630681233358.*i + 16083.45123273737) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10433310.88222013 - 625274.0027475406.*i))) + ((om1 - 1.0.*om + om2).*(- 0.0000000001143464823833147.*i - 0.00000002908519946984528))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(16083.45123273737 - 104.2630681233358.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (625274.0027475406.*i + 10433310.88222013))) + ((om1 - 1.0.*om + om2).*(- 0.0000000001582889038769531.*i - 0.0000000006499788787790134))./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(- 104.2630681233358.*i - 18580.67324419952) + om2.*(- 64.43133890242872.*i - 17956.36774133398) + om1.*om2 + om1.^2 + (755458.0170976262.*i + 11207692.7807475))));