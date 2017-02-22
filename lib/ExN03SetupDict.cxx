//
// File generated by rootcint at Thu Sep 10 15:04:49 2015

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME ExN03SetupDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "ExN03SetupDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void ExN03Setup_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_ExN03Setup(void *p = 0);
   static void *newArray_ExN03Setup(Long_t size, void *p);
   static void delete_ExN03Setup(void *p);
   static void deleteArray_ExN03Setup(void *p);
   static void destruct_ExN03Setup(void *p);
   static void streamer_ExN03Setup(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExN03Setup*)
   {
      ::ExN03Setup *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExN03Setup >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ExN03Setup", ::ExN03Setup::Class_Version(), "./ExN03Setup.hh", 16,
                  typeid(::ExN03Setup), DefineBehavior(ptr, ptr),
                  &::ExN03Setup::Dictionary, isa_proxy, 0,
                  sizeof(::ExN03Setup) );
      instance.SetNew(&new_ExN03Setup);
      instance.SetNewArray(&newArray_ExN03Setup);
      instance.SetDelete(&delete_ExN03Setup);
      instance.SetDeleteArray(&deleteArray_ExN03Setup);
      instance.SetDestructor(&destruct_ExN03Setup);
      instance.SetStreamerFunc(&streamer_ExN03Setup);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExN03Setup*)
   {
      return GenerateInitInstanceLocal((::ExN03Setup*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExN03Setup*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ExN03Setup::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ExN03Setup::Class_Name()
{
   return "ExN03Setup";
}

//______________________________________________________________________________
const char *ExN03Setup::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExN03Setup*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ExN03Setup::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExN03Setup*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void ExN03Setup::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExN03Setup*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *ExN03Setup::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gCINTMutex); if(!fgIsA) {fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExN03Setup*)0x0)->GetClass();} }
   return fgIsA;
}

//______________________________________________________________________________
void ExN03Setup::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExN03Setup.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> TargetRadius;
      R__b >> TargetLength;
      R__b >> ChamberInnerRadius;
      R__b >> ChamberThickness;
      R__b >> ChamberLength;
      R__b >> InnerRohacellThickness;
      R__b >> KaptonThickness;
      R__b >> OuterRohacellThickness;
      R__b >> TPCRadiusExt;
      R__b >> WindowThickness;
      R__b.CheckByteCount(R__s, R__c, ExN03Setup::IsA());
   } else {
      R__c = R__b.WriteVersion(ExN03Setup::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << TargetRadius;
      R__b << TargetLength;
      R__b << ChamberInnerRadius;
      R__b << ChamberThickness;
      R__b << ChamberLength;
      R__b << InnerRohacellThickness;
      R__b << KaptonThickness;
      R__b << OuterRohacellThickness;
      R__b << TPCRadiusExt;
      R__b << WindowThickness;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void ExN03Setup::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class ExN03Setup.
      TClass *R__cl = ::ExN03Setup::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "TargetRadius", &TargetRadius);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "TargetLength", &TargetLength);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "ChamberInnerRadius", &ChamberInnerRadius);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "ChamberThickness", &ChamberThickness);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "ChamberLength", &ChamberLength);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "InnerRohacellThickness", &InnerRohacellThickness);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "KaptonThickness", &KaptonThickness);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "OuterRohacellThickness", &OuterRohacellThickness);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "TPCRadiusExt", &TPCRadiusExt);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "WindowThickness", &WindowThickness);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExN03Setup(void *p) {
      return  p ? new(p) ::ExN03Setup : new ::ExN03Setup;
   }
   static void *newArray_ExN03Setup(Long_t nElements, void *p) {
      return p ? new(p) ::ExN03Setup[nElements] : new ::ExN03Setup[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExN03Setup(void *p) {
      delete ((::ExN03Setup*)p);
   }
   static void deleteArray_ExN03Setup(void *p) {
      delete [] ((::ExN03Setup*)p);
   }
   static void destruct_ExN03Setup(void *p) {
      typedef ::ExN03Setup current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_ExN03Setup(TBuffer &buf, void *obj) {
      ((::ExN03Setup*)obj)->::ExN03Setup::Streamer(buf);
   }
} // end of namespace ROOT for class ::ExN03Setup

/********************************************************
* ExN03SetupDict.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableExN03SetupDict();

extern "C" void G__set_cpp_environmentExN03SetupDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("ExN03Setup.hh");
  G__cpp_reset_tagtableExN03SetupDict();
}
#include <new>
extern "C" int G__cpp_dllrevExN03SetupDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* ExN03Setup */
static int G__ExN03SetupDict_184_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   ExN03Setup* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new ExN03Setup[n];
     } else {
       p = new((void*) gvp) ExN03Setup[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new ExN03Setup;
     } else {
       p = new((void*) gvp) ExN03Setup;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((ExN03Setup*) G__getstructoffset())->ReadConfigurationFile(*((string*) G__int(libp->para[0])));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ExN03Setup::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) ExN03Setup::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) ExN03Setup::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ExN03Setup::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((ExN03Setup*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) ExN03Setup::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ExN03Setup::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) ExN03Setup::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ExN03SetupDict_184_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ExN03Setup::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__ExN03SetupDict_184_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   ExN03Setup* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new ExN03Setup(*(ExN03Setup*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef ExN03Setup G__TExN03Setup;
static int G__ExN03SetupDict_184_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (ExN03Setup*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((ExN03Setup*) (soff+(sizeof(ExN03Setup)*i)))->~G__TExN03Setup();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (ExN03Setup*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((ExN03Setup*) (soff))->~G__TExN03Setup();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__ExN03SetupDict_184_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   ExN03Setup* dest = (ExN03Setup*) G__getstructoffset();
   *dest = *(ExN03Setup*) libp->para[0].ref;
   const ExN03Setup& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* ExN03Setup */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncExN03SetupDict {
 public:
  G__Sizep2memfuncExN03SetupDict(): p(&G__Sizep2memfuncExN03SetupDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncExN03SetupDict::*p)();
};

size_t G__get_sizep2memfuncExN03SetupDict()
{
  G__Sizep2memfuncExN03SetupDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceExN03SetupDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup))) {
     ExN03Setup *G__Lderived;
     G__Lderived=(ExN03Setup*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup),G__get_linked_tagnum(&G__ExN03SetupDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableExN03SetupDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__ExN03SetupDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__ExN03SetupDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__ExN03SetupDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__ExN03SetupDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__ExN03SetupDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__ExN03SetupDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__ExN03SetupDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__ExN03SetupDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__ExN03SetupDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__ExN03SetupDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* ExN03Setup */
static void G__setup_memvarExN03Setup(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup));
   { ExN03Setup *p; p=(ExN03Setup*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->TargetRadius)-(long)(p)),100,0,0,-1,-1,-1,1,"TargetRadius=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->TargetLength)-(long)(p)),100,0,0,-1,-1,-1,1,"TargetLength=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ChamberInnerRadius)-(long)(p)),100,0,0,-1,-1,-1,1,"ChamberInnerRadius=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ChamberThickness)-(long)(p)),100,0,0,-1,-1,-1,1,"ChamberThickness=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ChamberLength)-(long)(p)),100,0,0,-1,-1,-1,1,"ChamberLength=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->InnerRohacellThickness)-(long)(p)),100,0,0,-1,-1,-1,1,"InnerRohacellThickness=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->KaptonThickness)-(long)(p)),100,0,0,-1,-1,-1,1,"KaptonThickness=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->OuterRohacellThickness)-(long)(p)),100,0,0,-1,-1,-1,1,"OuterRohacellThickness=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->TPCRadiusExt)-(long)(p)),100,0,0,-1,-1,-1,1,"TPCRadiusExt=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->WindowThickness)-(long)(p)),100,0,0,-1,-1,-1,1,"WindowThickness=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__ExN03SetupDictLN_TClass),G__defined_typename("atomic_TClass_ptr"),-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarExN03SetupDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncExN03Setup(void) {
   /* ExN03Setup */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup));
   G__memfunc_setup("ExN03Setup",895,G__ExN03SetupDict_184_0_1, 105, G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("ReadConfigurationFile",2132,G__ExN03SetupDict_184_0_2, 121, -1, -1, 0, 1, 1, 1, 0, "u 'string' - 0 - -", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__ExN03SetupDict_184_0_3, 85, G__get_linked_tagnum(&G__ExN03SetupDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&ExN03Setup::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__ExN03SetupDict_184_0_4, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&ExN03Setup::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__ExN03SetupDict_184_0_5, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&ExN03Setup::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__ExN03SetupDict_184_0_6, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&ExN03Setup::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__ExN03SetupDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__ExN03SetupDict_184_0_10, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__ExN03SetupDict_184_0_11, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&ExN03Setup::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__ExN03SetupDict_184_0_12, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&ExN03Setup::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__ExN03SetupDict_184_0_13, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&ExN03Setup::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__ExN03SetupDict_184_0_14, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&ExN03Setup::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("ExN03Setup", 895, G__ExN03SetupDict_184_0_15, (int) ('i'), G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup), -1, 0, 1, 1, 1, 0, "u 'ExN03Setup' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~ExN03Setup", 1021, G__ExN03SetupDict_184_0_16, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__ExN03SetupDict_184_0_17, (int) ('u'), G__get_linked_tagnum(&G__ExN03SetupDictLN_ExN03Setup), -1, 1, 1, 1, 1, 0, "u 'ExN03Setup' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncExN03SetupDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalExN03SetupDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcExN03SetupDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__ExN03SetupDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__ExN03SetupDictLN_ExN03Setup = { "ExN03Setup" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableExN03SetupDict() {
  G__ExN03SetupDictLN_TClass.tagnum = -1 ;
  G__ExN03SetupDictLN_TBuffer.tagnum = -1 ;
  G__ExN03SetupDictLN_TMemberInspector.tagnum = -1 ;
  G__ExN03SetupDictLN_TObject.tagnum = -1 ;
  G__ExN03SetupDictLN_string.tagnum = -1 ;
  G__ExN03SetupDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__ExN03SetupDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__ExN03SetupDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__ExN03SetupDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__ExN03SetupDictLN_ExN03Setup.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableExN03SetupDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_string);
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__ExN03SetupDictLN_ExN03Setup),sizeof(ExN03Setup),-1,29952,(char*)NULL,G__setup_memvarExN03Setup,G__setup_memfuncExN03Setup);
}
extern "C" void G__cpp_setupExN03SetupDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupExN03SetupDict()");
  G__set_cpp_environmentExN03SetupDict();
  G__cpp_setup_tagtableExN03SetupDict();

  G__cpp_setup_inheritanceExN03SetupDict();

  G__cpp_setup_typetableExN03SetupDict();

  G__cpp_setup_memvarExN03SetupDict();

  G__cpp_setup_memfuncExN03SetupDict();
  G__cpp_setup_globalExN03SetupDict();
  G__cpp_setup_funcExN03SetupDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncExN03SetupDict();
  return;
}
class G__cpp_setup_initExN03SetupDict {
  public:
    G__cpp_setup_initExN03SetupDict() { G__add_setup_func("ExN03SetupDict",(G__incsetup)(&G__cpp_setupExN03SetupDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initExN03SetupDict() { G__remove_setup_func("ExN03SetupDict"); }
};
G__cpp_setup_initExN03SetupDict G__cpp_setup_initializerExN03SetupDict;
