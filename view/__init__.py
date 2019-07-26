from kivy.app import App
from kivy.lang import Builder
from kivy.uix.screenmanager import ScreenManager, Screen


Builder.load_string("""
<MenuScreen>:
    BoxLayout:
        canvas:
            Color:
                rgba: .7, .7, .75, 1
            Rectangle:
                pos: self.pos
                size: self.size
                 
        padding: 50
        spacing: 40
        Button:                                   
            text: 'Configurações'
            on_press: root.manager.current = 'config'
        Button:
            text: 'Sair'

<SettingsScreen>:
    FloatLayout:
        canvas:
            Color:
                rgba: .7, .7, .75, 1
            Rectangle:
                pos: self.pos
                size: self.size
                 
        padding: 50
        spacing: 40
        Label:
            text: "Experiment Configuration"
            pos_hint: {'center_x': 0.5, 'center_y': 0.9}
            font_size: 18

        Label:
            text: "Experiment Name"
            pos_hint: {'center_x':0.1, 'center_y':0.8}
            font_size: 14

        TextInput:
            id: name
            pos_hint: {'center_x':0.5, 'center_y':0.8}
            size_hint: 0.6, 0.07
            font_size: 14
            write_tab: False
            multiline: False

        Label:
            text: "Replicate Number"
            pos_hint: {'center_x':0.1, 'center_y':0.65}
            font_size: 14

        TextInput:
            id: replicate
            pos_hint: {'center_x':0.28, 'center_y':0.65}
            size_hint: 0.15, 0.07
            font_size: 14
            multiline: False
            write_tab: False
            input_type: 'number'
            input_filter: 'float'

        Label:
            text: "Group Number"
            pos_hint: {'center_x':0.55, 'center_y':0.65}
            font_size: 14

        TextInput:
            id: group
            pos_hint: {'center_x':0.7, 'center_y':0.65}
            size_hint: 0.15, 0.07
            font_size: 14
            multiline: False
            write_tab: False
            input_type: 'number'
            input_filter: 'float'
            
        Label:
            text: "Estoque Mínimo"
            pos_hint: {'center_x':0.1, 'center_y':0.35}
            font_size: 14

        TextInput:
            id: p_estoque_min
            pos_hint: {'center_x':0.5, 'center_y':0.35}
            size_hint: 0.6, 0.07
            font_size: 14
            write_tab: False
            multiline: False

        Button:
            text: "Salvar"
            pos_hint:{'center_x': 0.8, 'center_y': 0.15}
            size_hint: 0.2, 0.1
            on_press: app.salva_produto()
            font_size: 14
            
        Button:
            text: "Menu"
            pos_hint:{'center_x': 0.4, 'center_y': 0.15}
            size_hint: 0.2, 0.1
            on_press: root.manager.current = 'menu'
            font_size: 14
        
        
        
        
        
""")

# Declare both screens
# Button:
#             text: 'Voltar para o menu'
#             on_press: root.manager.current = 'menu'
class MenuScreen(Screen):
    pass

class SettingsScreen(Screen):
    pass

# Create the screen manager
sm = ScreenManager()
sm.add_widget(MenuScreen(name='menu'))
sm.add_widget(SettingsScreen(name='config'))

class TestApp(App):

    def build(self):
        return sm

if __name__ == '__main__':
    TestApp().run()