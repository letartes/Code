; ─────────────────────────────────────────────────────────────────────────────
; Inno Setup script — Protein Liability Analyzer  (Auto-HOS Edition)
; ─────────────────────────────────────────────────────────────────────────────
; Prerequisites:
;   1. Run build_windows.bat first to produce the .exe
;   2. Install Inno Setup 6 from https://jrsoftware.org/isinfo.php
;   3. Open this file in Inno Setup and press F9 to compile the installer
; ─────────────────────────────────────────────────────────────────────────────

#define AppName      "Protein Liability Analyzer"
#define AppVersion   "2.0"
#define AppPublisher "Letarte Scientific Consulting"
#define AppURL       "https://www.letartescientific.com"
#define AppExeName   "ProteinLiabilityAnalyzer.exe"

[Setup]
AppId={{A1B2C3D4-E5F6-7890-ABCD-EF1234567890}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
AppPublisherURL={#AppURL}
AppSupportURL={#AppURL}
AppUpdatesURL={#AppURL}
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
AllowNoIcons=yes
; Output location — one level above the build\ folder
OutputDir=..\dist
OutputBaseFilename=ProteinLiabilityAnalyzer_Setup_v{#AppVersion}
; Compression
Compression=lzma2/ultra64
SolidCompression=yes
; Require Windows 10 or later
MinVersion=10.0
; Warn if not running as admin
PrivilegesRequired=admin
PrivilegesRequiredOverridesAllowed=dialog
WizardStyle=modern
; Uncomment to add a custom icon:
; SetupIconFile=..\assets\icon.ico

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon";    Description: "{cm:CreateDesktopIcon}";    GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked
Name: "startmenuicon";  Description: "Create Start Menu shortcut"; GroupDescription: "{cm:AdditionalIcons}"; Flags: checkedonce

[Files]
; The main executable built by PyInstaller
Source: "..\dist\{#AppExeName}"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
; Start Menu shortcut
Name: "{group}\{#AppName}";         Filename: "{app}\{#AppExeName}"; Comment: "Protein sequence liability analysis with HOS prediction"
Name: "{group}\Uninstall {#AppName}"; Filename: "{uninstallexe}"
; Desktop shortcut (optional — controlled by task above)
Name: "{autodesktop}\{#AppName}";   Filename: "{app}\{#AppExeName}"; Tasks: desktopicon

[Run]
; Offer to launch the app after install
Filename: "{app}\{#AppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(AppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Clean up any HTML reports or logs written by the app into the install folder
Type: filesandordirs; Name: "{app}"
