$(function() {

//	-------------------------------------------------------------------------------------------------------
//		THINGS TO DO ON PAGE LOAD
//	-------------------------------------------------------------------------------------------------------

		// Copy global nav links into mobile menu
		
		if ( $('#globalLinks').length ) {
			$('#globalLinksMobile').html( $('#globalLinks').html() );
		} else {
			$('#globalLinksMobileHead').remove();
		}

		
		// Copy site nav links into mobile menu
		
		if ( $('#siteMenu').length ) {
			$('#siteMenuMobile').html( $('#siteMenu').html() );
		} else {
			//$('#siteMenuMobileHead').remove();
		}


		// Copy subnav (group nav) links into mobile menu
		
		if ( $('#subnav').length ) {
			$('#subnavMobile').html( $('#subnav').html() );
		} else {
			$('#subnavMobileHead').remove();
		}
		

		// Copy search to mobile search
	
		$('#siteMenuContainer').after('<div id="mobileSearch"></div>');
		$('#mobileSearch').html($('#gd-8').html());
		
		
		// Add megaMenu Class to Mega Menu items in site nav
		
		if( $('#siteMenu').length ) {
			$("#siteMenu > li").each(function(index) {
				$(this).attr('id','siteMenuItem'+(index+1));
				if ( $(this).find('.container_24').length > 0 ) {			
					$(this).children('ul').addClass("megaMenu");
				}
			});
		}


		// Add down caret to all top-level site menu <li> elements that have child drop-down lists
		
		var caret = "<span class='fa fa-caret-down menu-caret'></span>";
		$("#siteMenu > li").each(function(idx, li) {
			if( $(this).children('p').length > 0 ) {
				if ($(this).find('ul').length) {
					if( $(this).children('p').children('.fa').length === 0 ) {
						//alert( $(this).children('p').html() );
						$(this).children('p').append(caret);
					}
				}
			}
			if( $(this).children('a').length > 0 ) {
				if ($(this).find('ul').length) {
					if( $(this).children('a').children('.fa').length === 0 ) {
						//alert( $(this).children('a').html() );
						$(this).children('a').append(caret);
					}
				}
			}
		});
		
		$("#siteMenu").css("opacity", "1");


//	-------------------------------------------------------------------------------------------------------
//		GLOBAL NAV LINKS - ONCLICK FUNCTION
//	-------------------------------------------------------------------------------------------------------

		$('.gd-link').click(function() {
			console.log( "Garage door link clicked in header." );
			var linkIndex = $(this).parent('li').index();
			var GDid = "#gd-" + linkIndex;
			var GDdisplay = $(GDid).css('display');
			if (GDdisplay == "block") {
				$('.gd').slideUp('1000');
			} else {
				$.sidr('close','BNL-mobile-menu');
				$('.gd').slideUp('1000');
				$(GDid).slideDown('1000');
			}
		});
		
		$('.closeGD').click(function() {
				$('.gd').slideUp('1000');
		});


//	-------------------------------------------------------------------------------------------------------
//		MOBILE MENU (.sidr)
//	-------------------------------------------------------------------------------------------------------

		// Click hamburger icon to open mobile menu
		
		$('#hamburger').sidr({
			name:	'BNL-mobile-menu',
			displace: false
		});


		// Close Mobile Menu When User Clicks on Page Content
		
		$("html").on("click",function(e) {
			if ($(e.target).is('.sidr *, .sidr')) {
			
			} else {
				$.sidr('close','BNL-mobile-menu');
			}
		});


		// Close Mobile Menu Icon Click Function
		
		$(document).on("click", ".closeMobileMenuIcon", function(e) {
			$.sidr('close','BNL-mobile-menu');
		});
		
		$("#BNL-mobile-menu").on("swipeleft",function(){
			$.sidr('close', 'BNL-mobile-menu');
		});
		
		// jQuery Finger to Close Mobile Menu	
			
		$('body').on('flick', '#BNL-mobile-menu', function() {
			$.sidr('close', 'BNL-mobile-menu');
		});
		

//	-------------------------------------------------------------------------------------------------------
//		MOBILE SEARCH TRIGGER
//	-------------------------------------------------------------------------------------------------------

		$("#mobileSearchTrigger").click(function() {
			$('#mobileSearch').slideToggle();
		});			

			
//	-------------------------------------------------------------------------------------------------------
//		TOP OF PAGE TRIGGER
//	-------------------------------------------------------------------------------------------------------

		$(".topOfPageTrigger").click(function() {
			$.sidr('close','BNL-mobile-menu');
			$(".gd").hide();
			$('#mobileSearch').hide();
			$('html, body').animate({scrollTop: $("#header").offset().top}, 1000);
		});	
				
			
//	-------------------------------------------------------------------------------------------------------
//		DEPARTMENTS LINK (IN FOOTER)
//	-------------------------------------------------------------------------------------------------------

		$("#deptLink").click(function() {
			$.sidr('close','BNL-mobile-menu');
			$(".gd").hide();
			$('#mobileSearch').hide();
			$('html, body').animate({scrollTop: $("#header").offset().top}, 0);
			$('#gd-2').slideDown('1000');
		});		
		

//	-------------------------------------------------------------------------------------------------------
//		DISABLED BUTTONS (PREVENT CLICK)
//	-------------------------------------------------------------------------------------------------------
		
		$('body').on('click', '.button.disabled, #subnav .disabled, #subnav a.disabled', function(e) {
			e.preventDefault();
		});



//	WE MAY WANT TO MOVE THE FOLLOWING CODE INTO SEPARATE FILES SO WE ONLY INCLUDE IT ON PAGES THAT HAVE
//  THE ELEMENTS THAT REQUIRE IT

//	-------------------------------------------------------------------------------------------------------
//		TABS
//	-------------------------------------------------------------------------------------------------------

		// TABS WITH #tabs (ID, NOT CLASS)
		
		$( "#tabs" ).tabs();
		//REMEMBER SELECTED JQUERY TAB
		var tabCookieName = 'SelectedTab',
		$tabs = $('#tabs'),
		$lis = $tabs.find('ul').eq(0).find('li');
		
		$tabs.tabs({
			active: ( $.cookie( tabCookieName ) || 0 ),
			activate: function( e, ui ) {
				$.cookie( tabCookieName, $lis.index(ui.newTab) );
			}
		});

		// TABS WITH .tabs (CLASS, NOT ID)
		
		$( ".tabs" ).tabs();
		/*
		//REMEMBER SELECTED JQUERY TAB
		var tabCookieName = 'SelectedTab2',
		$tabs = $('.tabs'),
		$lis = $tabs.find('ul').eq(0).find('li');
		
		$tabs.tabs({
			active: ( $.cookie( tabCookieName ) || 0 ),
			activate: function( e, ui ) {
				$.cookie( tabCookieName, $lis.index(ui.newTab) );
			}
		});
		*/


//	-------------------------------------------------------------------------------------------------------
//		VERTICAL TABS - REMEMBER SELECTED JQUERY TAB
//	-------------------------------------------------------------------------------------------------------
		
		var vtabCookieName = 'SelectedTab',
			$vtabs = $('#vertical-tabs'),
			$vlis = $vtabs.find('ul').eq(0).find('li');
		
		$vtabs.tabs({
			active: ( $.cookie( vtabCookieName ) || 0 ),
			activate: function( e, ui ) {
				$.cookie( vtabCookieName, $vlis.index(ui.newTab) );
			}
		});



});