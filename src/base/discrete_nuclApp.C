#include "discrete_nuclApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
discrete_nuclApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

discrete_nuclApp::discrete_nuclApp(InputParameters parameters) : MooseApp(parameters)
{
  discrete_nuclApp::registerAll(_factory, _action_factory, _syntax);
}

discrete_nuclApp::~discrete_nuclApp() {}

void
discrete_nuclApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"discrete_nuclApp"});
  Registry::registerActionsTo(af, {"discrete_nuclApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
discrete_nuclApp::registerApps()
{
  registerApp(discrete_nuclApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
discrete_nuclApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  discrete_nuclApp::registerAll(f, af, s);
}
extern "C" void
discrete_nuclApp__registerApps()
{
  discrete_nuclApp::registerApps();
}
